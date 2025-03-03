/*
 *  VISDEM 1.2 
 *
 *  Copyright 2025 Gunnar F. Schroeder
 *
 *  VISDEM is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  VISDEM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  If not, see <http://www.gnu.org/licenses/>.
 * 
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <fftw3.h>


/* -------------------------
   Data Structures
   ------------------------- */

// A minimal structure to hold a 3D density map.
typedef struct {
    int nx, ny, nz;         // Grid dimensions (voxels)
    double apix_x, apix_y, apix_z;  // Pixel sizes (Å)
    double origin_x, origin_y, origin_z; // Origin coordinates (Å)
    int n;                  // Total number of voxels = nx*ny*nz
    double *data;           // Density data (stored as doubles)
} DensityMap;


/* Global constant for Peng kernel */
double PengGFactor = 10.0;

/*** Electron Scattering ****/
/*  Taken from Table 3 in 
  L.-M. PENG,G. REN, S. L. DUDAREV AND M. J. WHELAN
  Acta Cryst. (1996). A52, 257-276

   Coefficient arrays for the Peng kernel.
   Each row corresponds to an atom type. For example, H=1, C=2, N=3, O=4, S=5, etc.
*/
double PengScatterFactor_A[][5] = {
    { 0.0489, 0.2091, 0.7537, 1.1420, 0.3555 }, // ANY
    { 0.0088, 0.0449, 0.1481, 0.2356, 0.0914 }, // H
    { 0.0489, 0.2091, 0.7537, 1.1420, 0.3555 }, // C
    { 0.0267, 0.1328, 0.5301, 1.1020, 0.4215 }, // N
    { 0.0365, 0.1729, 0.5805, 0.8814, 0.3121 }, // O
    { 0.0915, 0.4312, 1.0847, 2.4671, 1.0852 }, // S
    { 0.1005, 0.4615, 1.0663, 2.5854, 1.2725 }, // P
    { 0.5229, 2.2874, 4.7243, 5.0807, 5.6389 }, // BA
    { 0.1260, 0.6442, 0.8893, 1.8197, 1.2988 }, // NA
    { 0.0799, 0.3891, 1.0037, 2.3332, 1.0507 }, // CL
    { 0.2149, 0.8703, 2.4999, 2.3591, 3.0318 }, // K
    { 0.1130, 0.5575, 0.9046, 2.1580, 1.4735 }  // MG
};

double PengScatterFactor_B[][5] = {
    { 0.1140, 1.0825, 5.4281, 17.8811, 51.1341 }, // ANY
    { 0.1152, 1.0867, 4.9755, 16.5591, 43.2743 }, // H
    { 0.1140, 1.0825, 5.4281, 17.8811, 51.1341 }, // C
    { 0.0541, 0.5165, 2.8207, 10.6297, 34.3764 }, // N
    { 0.0652, 0.6184, 2.9449, 9.6298, 28.2194 },  // O
    { 0.0838, 0.7788, 4.3462, 15.5846, 44.6365 }, // S
    { 0.0977, 0.9084, 4.9654, 18.5471, 54.3648 }, // P
    { 0.1434, 1.6019, 9.4511, 42.7685, 148.4969 }, // BA
    { 0.1684, 1.7150, 8.8386, 50.8265, 147.2073 }, // NA
    { 0.0694, 0.6443, 3.5351, 12.5058, 35.8633 },  // CL
    { 0.1660, 1.6906, 8.7447, 46.7825, 165.6923 }, // K
    { 0.1356, 1.3579, 6.9255, 32.3165, 92.1138 }   // MG
};

/*
 * PengKernel computes the density kernel for a given squared distance 'r',
 * an atom type (index into the coefficient arrays), and sigma.
 * It returns a sum of five Gaussian-like components.
 */
static double PengKernel(double r, int type, double sigma)
{
    double v = 0.0, h1, h2, h4, h0;
    double sf[5], sf_a[5];

    h0 = -9.8696 * r;  // constant factor times r

    for (int comp = 0; comp < 5; comp++) {
        sf[comp]   = PengScatterFactor_B[type][comp];
        sf_a[comp] = PengScatterFactor_A[type][comp];
    }
    for (int comp = 0; comp < 5; comp++) {
        h1 = sf[comp];
        h2 = 1.0 / (h1 + PengGFactor);
        h4 = pow(h2, 1.5);  // equivalent to sqrt(h2^3)
        v += sf_a[comp] * h4 * exp(h0 * h2);
    }
    return v;
}

/*
 * renderPeng renders a density map from an atomic model using the Peng kernel.
 *
 * Parameters:
 *   result: Pointer to a DensityMap structure (the output density map).
 *   width:  The spatial width (in Å) over which each atom contributes.
 *   sigma:  The sigma parameter for the kernel.
 *   n:      The number of atoms.
 *   type:   An integer array (length n) with atom type codes (e.g., H=1, C=2, etc.).
 *   coords: A double array (length 3*n) with the x,y,z coordinates for each atom.
 *   factor: A double array (length n) with a mass or scattering factor per atom.
 *
 * This function loops over all atoms, computes the grid region around each atom, and
 * adds the density contributions (via PengKernel) to the density map.
 */
DensityMap * renderPeng(DensityMap *result, double width, double sigma, size_t n, int *type, double *coords, double *factor)
//void renderPeng(DensityMap *result, double width, double sigma, size_t n, int *type, double *coords, double *factor)
{
    size_t atom;
    size_t i, j, k;
    size_t idx, idxx;
    size_t dim2 = result->nx;                 // number of grid points in x
    size_t dim3 = result->nx * result->ny;      // number of grid points per slice
    // Determine the number of grid points (walk) to cover the kernel width.
    size_t walk_x = (size_t)(0.5 * width / result->apix_x) + 1;
    size_t walk_y = (size_t)(0.5 * width / result->apix_y) + 1;
    size_t walk_z = (size_t)(0.5 * width / result->apix_z) + 1;

    // Loop over each atom.
    for (atom = 0; atom < n; atom++) {
        double mass = factor[atom];
        int type_code = type[atom];

        // Get atom position relative to map origin.
        double pos_x = coords[3 * atom    ] - result->origin_x;
        double pos_y = coords[3 * atom + 1] - result->origin_y;
        double pos_z = coords[3 * atom + 2] - result->origin_z;

        // Determine the grid index corresponding to the atom position.
        size_t gridid_x = (size_t)(pos_x / result->apix_x);
        size_t gridid_y = (size_t)(pos_y / result->apix_y);
        size_t gridid_z = (size_t)(pos_z / result->apix_z);

        // Set the extent of the kernel region.
        size_t gridend_x = gridid_x + walk_x;
        size_t gridend_y = gridid_y + walk_y;
        size_t gridend_z = gridid_z + walk_z;
        if (gridid_x > walk_x)
            gridid_x = gridid_x - walk_x + 1;
        else
            gridid_x = 0;
        if (gridid_y > walk_y)
            gridid_y = gridid_y - walk_y + 1;
        else
            gridid_y = 0;
        if (gridid_z > walk_z)
            gridid_z = gridid_z - walk_z + 1;
        else
            gridid_z = 0;

        // Loop over the small cube around the atom.
        for (k = gridid_z; k < gridend_z; k++) {
            if (k < (size_t)result->nz) {
                idxx = k * dim3;
                for (j = gridid_y; j < gridend_y; j++) {
                    if (j < (size_t)result->ny) {
                        for (i = gridid_x; i < gridend_x; i++) {
                            if (i < (size_t)result->nx) {
                                idx = i + idxx + j * dim2;
                                // Compute the squared distance from the grid point to the atom.
                                double rx = result->apix_x * i - pos_x;
                                double ry = result->apix_y * j - pos_y;
                                double rz = result->apix_z * k - pos_z;
                                double len2 = rx * rx + ry * ry + rz * rz;
                                // Evaluate the Peng kernel at this distance.
                                double dens = PengKernel(len2, type_code, sigma * sigma);
                                dens *= mass;
                                result->data[idx] += dens;
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}



/* -------------------------
   MRC I/O Functions
   ------------------------- */

/*
 * readMRC reads a cryo-EM MRC file.
 * It reads a 1024-byte header into the provided header buffer,
 * then extracts the grid dimensions, cell dimensions (to derive pixel size),
 * and origin from standard offsets.
 * It assumes that the image data are stored as 32-bit floats (mode 2).
 */
int readMRC(const char *filename, DensityMap *map, char header[1024]) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        return 1;
    }
    if (fread(header, 1, 1024, fp) != 1024) {
        fprintf(stderr, "Error reading header from %s\n", filename);
        fclose(fp);
        return 1;
    }
    int nx = *((int *) &header[0]);
    int ny = *((int *) &header[4]);
    int nz = *((int *) &header[8]);
    int mode = *((int *) &header[12]);
    if (mode != 2)
        fprintf(stderr, "Warning: mode != 2; only float32 (mode 2) is supported\n");
    map->nx = nx; map->ny = ny; map->nz = nz;
    map->n = nx * ny * nz;
    // Cell dimensions (in Å) are stored at offsets 40, 44, 48.
    float xlen = *((float *) &header[40]);
    float ylen = *((float *) &header[44]);
    float zlen = *((float *) &header[48]);
    map->apix_x = (nx > 0) ? xlen / nx : 1.0;
    map->apix_y = (ny > 0) ? ylen / ny : 1.0;
    map->apix_z = (nz > 0) ? zlen / nz : 1.0;
    // Origin stored at offset 196.
    float xorigin = *((float *) &header[196]);
    float yorigin = *((float *) &header[200]);
    float zorigin = *((float *) &header[204]);
    map->origin_x = xorigin;
    map->origin_y = yorigin;
    map->origin_z = zorigin;
    map->data = (double *) malloc(sizeof(double) * map->n);
    if (!map->data) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(fp);
        return 1;
    }
    float *temp = (float *) malloc(sizeof(float) * map->n);
    if (!temp) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(fp);
        return 1;
    }
    if (fread(temp, sizeof(float), map->n, fp) != (size_t) map->n) {
        fprintf(stderr, "Error reading map data from %s\n", filename);
        free(temp);
        fclose(fp);
        return 1;
    }
    for (int i = 0; i < map->n; i++) {
        map->data[i] = temp[i];
    }
    free(temp);
    fclose(fp);
    return 0;
}

/*
 * writeMRC writes the density map to an MRC file.
 * Before writing, it updates key header fields (dimensions, mode, cell lengths, origin)
 * in the provided 1024-byte header buffer.
 */
int writeMRC(const char *filename, DensityMap *map, char header[1024]) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Cannot open file %s for writing\n", filename);
        return 1;
    }
    int nx = map->nx, ny = map->ny, nz = map->nz;
    *((int *) &header[0]) = nx;
    *((int *) &header[4]) = ny;
    *((int *) &header[8]) = nz;
    *((int *) &header[12]) = 2; // mode 2 (float32)
    float xlen = map->apix_x * nx;
    float ylen = map->apix_y * ny;
    float zlen = map->apix_z * nz;
    *((float *) &header[40]) = xlen;
    *((float *) &header[44]) = ylen;
    *((float *) &header[48]) = zlen;
    *((float *) &header[196]) = (float) map->origin_x;
    *((float *) &header[200]) = (float) map->origin_y;
    *((float *) &header[204]) = (float) map->origin_z;
    if (fwrite(header, 1, 1024, fp) != 1024) {
        fprintf(stderr, "Error writing header to %s\n", filename);
        fclose(fp);
        return 1;
    }
    float *temp = (float *) malloc(sizeof(float) * map->n);
    if (!temp) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(fp);
        return 1;
    }
    for (int i = 0; i < map->n; i++) {
        temp[i] = (float) map->data[i];
    }
    if (fwrite(temp, sizeof(float), map->n, fp) != (size_t) map->n) {
        fprintf(stderr, "Error writing map data to %s\n", filename);
        free(temp);
        fclose(fp);
        return 1;
    }
    free(temp);
    fclose(fp);
    return 0;
}

/* -------------------------
   FFTW-Based Smoothing Filter
   ------------------------- */

/*
 * smoothLowpassFilter_many applies a smooth (Butterworth) low-pass filter
 * to the density map using FFTW’s “many” interface.
 */
void smoothLowpassFilter_many(DensityMap *map, double cutoff_resolution) {
    int nx = map->nx, ny = map->ny, nz = map->nz;
    int n = map->n;
    int rank = 3;
    int dims[3] = { nz, ny, nx };
    int howmany = 1;
    int inembed[3] = { nz, ny, nx };
    int onembed[3] = { nz, ny, nx/2 + 1 };
    int istride = 1, ostride = 1;
    int idist = n;
    int odist = nz * ny * (nx/2 + 1);
    
    fftw_complex *fft_data = fftw_malloc(sizeof(fftw_complex) * odist);
    if (!fft_data) {
        fprintf(stderr, "Error allocating FFT data\n");
        exit(1);
    }
    fftw_plan plan_forward = fftw_plan_many_dft_r2c(rank, dims, howmany,
                                    map->data, inembed, istride, idist,
                                    fft_data, onembed, ostride, odist,
                                    FFTW_ESTIMATE);
    fftw_execute(plan_forward);
    fftw_destroy_plan(plan_forward);
    
    double Lx = nx * map->apix_x;
    double Ly = ny * map->apix_y;
    double Lz = nz * map->apix_z;
    double f_cut = 1.0 / cutoff_resolution;
    int order = 4;
    for (int z = 0; z < nz; z++) {
        int kz = (z <= nz/2) ? z : z - nz;
        double freq_z = (double) kz / Lz;
        for (int y = 0; y < ny; y++) {
            int ky = (y <= ny/2) ? y : y - ny;
            double freq_y = (double) ky / Ly;
            for (int x = 0; x < (nx/2 + 1); x++) {
                double freq_x = (double) x / Lx;
                double freq_mag = sqrt(freq_x*freq_x + freq_y*freq_y + freq_z*freq_z);
                double factor = 1.0 / (1.0 + pow(freq_mag / f_cut, 2 * order));
                int index = z * ny * (nx/2 + 1) + y * (nx/2 + 1) + x;
                fft_data[index][0] *= factor;
                fft_data[index][1] *= factor;
            }
        }
    }
    fftw_plan plan_inverse = fftw_plan_many_dft_c2r(rank, dims, howmany,
                                    fft_data, onembed, ostride, odist,
                                    map->data, inembed, istride, idist,
                                    FFTW_ESTIMATE);
    fftw_execute(plan_inverse);
    fftw_destroy_plan(plan_inverse);
    fftw_free(fft_data);
    for (int i = 0; i < n; i++) {
        map->data[i] /= n;
    }
}

/* -------------------------
   Utility Functions
   ------------------------- */

// gaussRand generates a Gaussian random number (Box-Muller).
double gaussRand() {
    double r1, r2, q, p;
    static double r2_cached = 0.0;
    static int have_cached = 0;
    if (have_cached) {
        have_cached = 0;
        return r2_cached;
    }
    do {
        r1 = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
        r2 = 2.0 * ((double) rand() / RAND_MAX) - 1.0;
        q = r1 * r1 + r2 * r2;
    } while(q == 0.0 || q >= 1.0);
    p = sqrt(-2.0 * log(q) / q);
    r2_cached = r2 * p;
    have_cached = 1;
    return r1 * p;
}

// cmpfunc is used by qsort to sort double pointers in ascending order.
int cmpfunc(const void *a, const void *b) {
    double da = **(double **)a;
    double db = **(double **)b;
    return (da < db) ? -1 : (da > db);
}


/* 
 * matchStructureFactor adjusts the Fourier amplitudes of inmap so that the
 * radial (i.e. spherically averaged) amplitude spectrum of inmap matches that of beadMap.
 * Both maps must have the same dimensions.
 */
void matchStructureFactor(DensityMap *inmap, DensityMap *beadMap, int nbins) {
    int nx = inmap->nx, ny = inmap->ny, nz = inmap->nz;
    int nxh = nx/2 + 1;
    int n = inmap->n;
    
    // Allocate arrays for Fourier transforms.
    fftw_complex *F_in = fftw_malloc(sizeof(fftw_complex) * nz * ny * nxh);
    fftw_complex *F_bead = fftw_malloc(sizeof(fftw_complex) * nz * ny * nxh);
    if (!F_in || !F_bead) {
        fprintf(stderr, "Memory allocation error in matchStructureFactor\n");
        exit(1);
    }
    
    // Compute forward FFTs.
    fftw_plan plan_in = fftw_plan_dft_r2c_3d(nz, ny, nx, inmap->data, F_in, FFTW_ESTIMATE);
    fftw_plan plan_bead = fftw_plan_dft_r2c_3d(nz, ny, nx, beadMap->data, F_bead, FFTW_ESTIMATE);
    fftw_execute(plan_in);
    fftw_execute(plan_bead);
    fftw_destroy_plan(plan_in);
    fftw_destroy_plan(plan_bead);
    
    // Physical dimensions.
    double Lx = nx * inmap->apix_x;
    double Ly = ny * inmap->apix_y;
    double Lz = nz * inmap->apix_z;
    // Maximum frequency (approximate, from Nyquist frequencies)
    double max_freq = sqrt(pow((nx/2.0)/Lx,2) + pow((ny/2.0)/Ly,2) + pow((nz/2.0)/Lz,2));
    double bin_size = max_freq / nbins;
    
    // Allocate arrays to accumulate amplitudes and counts.
    double *sum_in   = calloc(nbins, sizeof(double));
    double *sum_bead = calloc(nbins, sizeof(double));
    int    *count_in   = calloc(nbins, sizeof(int));
    int    *count_bead = calloc(nbins, sizeof(int));
    if (!sum_in || !sum_bead || !count_in || !count_bead) {
        fprintf(stderr, "Memory allocation error in matchStructureFactor\n");
        exit(1);
    }
    
    // Loop over Fourier coefficients (for both inmap and beadMap).
    for (int z = 0; z < nz; z++) {
        int kz = (z <= nz/2) ? z : z - nz; // account for negative frequencies
        double fz = (double) kz / Lz;
        for (int y = 0; y < ny; y++) {
            int ky = (y <= ny/2) ? y : y - ny;
            double fy = (double) ky / Ly;
            for (int x = 0; x < nxh; x++) {
                double fx = (double) x / Lx;
                double freq = sqrt(fx*fx + fy*fy + fz*fz);
                int bin = (int)(freq / bin_size);
                if (bin >= nbins) bin = nbins - 1;
                int index = z * ny * nxh + y * nxh + x;
                double amp_in   = sqrt(F_in[index][0]*F_in[index][0] + F_in[index][1]*F_in[index][1]);
                double amp_bead = sqrt(F_bead[index][0]*F_bead[index][0] + F_bead[index][1]*F_bead[index][1]);
                sum_in[bin]   += amp_in;
                sum_bead[bin] += amp_bead;
                count_in[bin]++;
                count_bead[bin]++;
            }
        }
    }
    
    // Compute the radially averaged amplitudes.
    double *A_in = malloc(nbins * sizeof(double));
    double *A_bead = malloc(nbins * sizeof(double));
    for (int i = 0; i < nbins; i++) {
        A_in[i] = (count_in[i] > 0) ? (sum_in[i] / count_in[i]) : 0.0;
        A_bead[i] = (count_bead[i] > 0) ? (sum_bead[i] / count_bead[i]) : 0.0;
    }
    
    // Now, for each Fourier coefficient in inmap, compute a scale factor so that
    // its amplitude is adjusted in proportion to the ratio A_bead / A_in for its radial bin.
    for (int z = 0; z < nz; z++) {
        int kz = (z <= nz/2) ? z : z - nz;
        double fz = (double) kz / Lz;
        for (int y = 0; y < ny; y++) {
            int ky = (y <= ny/2) ? y : y - ny;
            double fy = (double) ky / Ly;
            for (int x = 0; x < nxh; x++) {
                double fx = (double) x / Lx;
                double freq = sqrt(fx*fx + fy*fy + fz*fz);
                int bin = (int)(freq / bin_size);
                if (bin >= nbins) bin = nbins - 1;
                double factor = 1.0;
                if (A_in[bin] > 1e-6)
                    factor = A_bead[bin] / A_in[bin];
                int index = z * ny * nxh + y * nxh + x;
                F_in[index][0] *= factor;
                F_in[index][1] *= factor;
            }
        }
    }
    
    // Compute the inverse FFT to update inmap->data.
    fftw_plan plan_inv = fftw_plan_dft_c2r_3d(nz, ny, nx, F_in, inmap->data, FFTW_ESTIMATE);
    fftw_execute(plan_inv);
    fftw_destroy_plan(plan_inv);
    // Normalize inverse FFT.
    for (int i = 0; i < n; i++) {
        inmap->data[i] /= n;
    }
    
    fftw_free(F_in);
    fftw_free(F_bead);
    free(sum_in);
    free(sum_bead);
    free(count_in);
    free(count_bead);
    free(A_in);
    free(A_bead);
}


// A helper struct to store a voxel’s value together with its index.
typedef struct {
    int index;
    double value;
} IndexValue;

// Comparator for qsort: sorts in ascending order based on the value.
int cmpIndexValue(const void *a, const void *b) {
    const IndexValue *ia = (const IndexValue *) a;
    const IndexValue *ib = (const IndexValue *) b;
    if (ia->value < ib->value)
        return -1;
    else if (ia->value > ib->value)
        return 1;
    else
        return 0;
}

/*
 * matchHistogram replaces the voxel values in target so that when sorted,
 * they are exactly identical to the voxel values in source.
 *
 * This is done by:
 *  1. Creating an array of IndexValue for target and source.
 *  2. Sorting both arrays in ascending order by value.
 *  3. For each rank i, assigning:
 *         target[ sorted_target[i].index ] = sorted_source[i].value;
 */
void matchHistogram(DensityMap *target, DensityMap *source) {
    int n = target->n;
    IndexValue *tArr = malloc(n * sizeof(IndexValue));
    IndexValue *sArr = malloc(n * sizeof(IndexValue));
    if (!tArr || !sArr) {
        fprintf(stderr, "Memory allocation error in matchHistogram\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        tArr[i].index = i;
        tArr[i].value = target->data[i];
        sArr[i].index = i;
        sArr[i].value = source->data[i];
    }
    qsort(tArr, n, sizeof(IndexValue), cmpIndexValue);
    qsort(sArr, n, sizeof(IndexValue), cmpIndexValue);
    // Now reassign target values: for each rank i, copy the source value.
    for (int i = 0; i < n; i++) {
        target->data[tArr[i].index] = sArr[i].value;
    }
    free(tArr);
    free(sArr);
}

// Write the bead model (array of coordinates) in PDB format.
void writePDB(const char *filename, double *xcoo, int natoms) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s for writing\n", filename);
        return;
    }
    for (int i = 0; i < natoms; i++) {
        fprintf(fp,
            "ATOM  %5d  CA  BEA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
            i + 1, i + 1, xcoo[i * 3 + 0], xcoo[i * 3 + 1], xcoo[i * 3 + 2]);
    }
    fprintf(fp, "END\n");
    fclose(fp);
}



/* -------------------------
   Main Application
   ------------------------- */

int main(int argc, char *argv[]) {

    char *inmap_fname = NULL;
    char *outmap_fname = NULL;
    double mass = -1.0;              // in kDa
    double resolution = 2.0;        // resolution cutoff in Å
    double threshold = -1.0;        // optional threshold
    double rise = 4.75, twist = -0.6; // optional helical parameters
    int use_threshold = 0, use_symmetry = 0;
    
    int c;
    int option_index = 0;
    
    static struct option long_options[] = {
        {"input",      required_argument, 0, 'i'},
        {"output",     required_argument, 0, 'o'},
        {"mass",       required_argument, 0, 'm'},
        {"resolution", required_argument, 0, 'r'},
        {"threshold",  required_argument, 0, 't'},
        {"rise",       required_argument, 0,  0 },
        {"twist",      required_argument, 0,  0 },
        {0, 0, 0, 0}
    };

    printf("VISDEM 1.2 (2025) by Gunnar F. Schroeder, Michaela Spiegel, and Amudha Duraisamy\n");


    
    while ((c = getopt_long(argc, argv, "i:o:m:r:t:", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i': 
                inmap_fname = optarg; 
                break;
            case 'o': 
                outmap_fname = optarg; 
                break;
            case 'm': 
                mass = atof(optarg); 
                break;
            case 'r': 
                resolution = atof(optarg); 
                break;
            case 't': 
                threshold = atof(optarg);
                use_threshold = 1;
                break;
            case 0:
                if (strcmp(long_options[option_index].name, "rise") == 0) {
                    rise = atof(optarg);
                    use_symmetry = 1;
                } else if (strcmp(long_options[option_index].name, "twist") == 0) {
                    twist = atof(optarg);
                    use_symmetry = 1;
                }
                break;
            default:
                fprintf(stderr, "Usage: %s --input <file> --output <file> --mass <mass> --resolution <cutoff> [--threshold <value>] [--rise <value>] [--twist <value>]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    // Check required parameters.
    if (!inmap_fname || !outmap_fname || (mass <= 0.0 && !use_threshold)) {
        fprintf(stderr, "Missing required parameters.\n");
        fprintf(stderr, "Usage: %s --input <file> --output <file> --mass <mass> --resolution <cutoff> [--threshold <value>] [--rise <value> --twist <value>]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    
    // Debug output.
    printf("Input map: %s\n", inmap_fname);
    printf("Output map: %s\n", outmap_fname);
    //printf("Mass (kDa): %.2lf\n", mass);
    //printf("Resolution cutoff (Å): %.2lf\n", resolution);

    /*
    if (use_threshold)
        printf("Threshold: %lf\n", threshold);
    if (use_symmetry) {
        printf("Rise: %.3lf\n", rise);
        printf("Twist: %.3lf\n", twist);
    }
    */
    
    mass *= 1000.0; // Convert mass from kDa to Da.



    DensityMap inmap;
    char header[1024];
    if (readMRC(inmap_fname, &inmap, header)) {
        fprintf(stderr, "ERROR reading %s\n", inmap_fname);
        exit(1);
    }
    
    // Apply a smooth low-pass filter to the input map.
    // TODO: decide maybe apply filter also to input ?
    //smoothLowpassFilter_many(&inmap, resolution);
    
    // Compute basic statistics.
    double min_val = 1e9, max_val = -1e9, sum_val = 0.0;
    for (int i = 0; i < inmap.n; i++) {
        double v = inmap.data[i];
        sum_val += v;
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
    }
    double mean_val = sum_val / inmap.n;
    fprintf(stderr, "Map: min = %lf, max = %lf, mean = %lf\n", min_val, max_val, mean_val);

    double voxel_volume = inmap.apix_x * inmap.apix_y * inmap.apix_z;
    double avg_mass = 1.0*0.499 + 12.0*0.314 + 14.0*0.087 + 16.0*0.098 + 32.0*0.002;
    int natoms = 0;
    int nvoxels = 0;
    double volume = 0.0;
    
    // If no threshold is provided, determine one from the density histogram.
    if (!use_threshold && mass > 0) {
        /* Case 1: Only mass is provided.
        Compute the total protein volume from mass, then compute the number of voxels 
        that should be above threshold, and derive threshold from the sorted density values.
        */
        // Compute number of voxels corresponding to the protein volume.
        volume = mass * 1.19; // Å^3 per Da
        nvoxels = (int)(volume / voxel_volume);
        double **values = (double **) malloc(inmap.n * sizeof(double *));
        if (!values) { fprintf(stderr, "Memory allocation error\n"); exit(1); }
        for (int i = 0; i < inmap.n; i++) {
            values[i] = &inmap.data[i];
        }
        qsort(values, inmap.n, sizeof(double *), cmpfunc);
        int idx = inmap.n - nvoxels;
        if (idx < 0) idx = 0;
        threshold = *(values[idx]);
        free(values);
        natoms = (int)(mass / avg_mass);
    } else if (mass <=0 &&  use_threshold) {
        /* Case 2: Only threshold is provided.
        Count the number of voxels above threshold, compute the volume above threshold,
        and then estimate natoms from an assumed average atomic volume (e.g., 15 Å³/atom).
        */
        for (int i = 0; i < inmap.n; i++) {
           if (inmap.data[i] >= threshold)
               nvoxels++;
        }
        volume = nvoxels * voxel_volume;
        mass = volume / 1.19;
        natoms = (int)(mass / avg_mass);

    } else if (mass > 0 && use_threshold) {
        /* Case 3: Both mass and threshold are provided.
        Use the provided threshold and compute natoms from mass.
        */
        natoms = (int)(mass / avg_mass);
    }

    //double volume_above_threshold = nvoxels * voxel_volume;
    //natoms = (int)(volume_above_threshold / 15.0);

    printf("Mass (kDa): %.2lf\n", mass/1000.0);
    printf("Volume (A^3): %.2lf\n", volume);
    printf("Using density threshold = %lf\n", threshold);
    

    // --- Compute Volume and Atom Counts ---
    /* Given:
       volume = mass * 1.19;  // in Å^3 (using 1.19 Å^3 per Da)
       and average atomic fractions:
         H = 0.499
         C = 0.314
         N = 0.087
         O = 0.098
         S = 0.002
       We compute the effective average mass per atom:
    */
    int natoms_c = (int)(natoms * 0.314);
    int natoms_n = (int)(natoms * 0.087);
    int natoms_o = (int)(natoms * 0.098);
    int natoms_s = (int)(natoms * 0.002);
    int natoms_h = natoms - natoms_c - natoms_n - natoms_o - natoms_s;

    fprintf(stderr, "Total atoms = %d, C = %d, N = %d, O = %d, S = %d, H = %d\n",
        natoms, natoms_c, natoms_n, natoms_o, natoms_s, natoms_h);

    // Allocate and fill an array of atom types (using codes: H=1, C=2, N=3, O=4, S=5).
    int *atom_type = (int *) malloc(natoms * sizeof(int));
    if (!atom_type) {
        fprintf(stderr, "Memory allocation error for atom_type\n");
        exit(1);
    }
    //int i, j, pos = 0;
    int i, pos = 0;
    for (i = 0; i < natoms_c; i++) {
        atom_type[pos++] = 2;  // Carbon
    }
    for (i = 0; i < natoms_n; i++) {
        atom_type[pos++] = 3;  // Nitrogen
    }
    for (i = 0; i < natoms_o; i++) {
        atom_type[pos++] = 4;  // Oxygen
    }
    for (i = 0; i < natoms_s; i++) {
        atom_type[pos++] = 5;  // Sulfur
    }
    for (i = 0; i < natoms_h; i++) {
        atom_type[pos++] = 1;  // Hydrogen
    }
    // Shuffle the atom_type array (Fisher-Yates shuffle).
    for (i = natoms - 1; i > 0; i--) {
        int k = rand() % (i + 1);
        int tmp = atom_type[i];
        atom_type[i] = atom_type[k];
        atom_type[k] = tmp;
    }
    

    // Allocate array for bead coordinates.
    double *xcoo = (double *) calloc(3 * natoms, sizeof(double));
    if (!xcoo) { fprintf(stderr, "Memory allocation error\n"); exit(1); }


    //
    // Build bead model: randomly choose voxels above threshold.
    //

    // Compute map boundaries.
    double x_min = inmap.origin_x;
    double x_max = inmap.origin_x + inmap.nx * inmap.apix_x;
    double y_min = inmap.origin_y;
    double y_max = inmap.origin_y + inmap.ny * inmap.apix_y;
    double z_min = inmap.origin_z;
    double z_max = inmap.origin_z + inmap.nz * inmap.apix_z;
    
    // Compute center of the map for rotation (assuming helix axis along z).
    double center_x = inmap.origin_x + (inmap.nx * inmap.apix_x) / 2.0;
    double center_y = inmap.origin_y + (inmap.ny * inmap.apix_y) / 2.0;


   double twist_rad = twist * (M_PI / 180.0);
   int j_max = 20; // Maximum number of symmetry steps per ASU bead
   
   int bead_count = 0;
   while (bead_count < natoms) {
       // Select a candidate ASU bead from a random voxel above threshold.
       int valid = 0, ind1, ind2, ind3, idx;
       double asu_x, asu_y, asu_z;
       while (!valid) {
           ind1 = rand() % inmap.nx;
           ind2 = rand() % inmap.ny;
           ind3 = rand() % inmap.nz;
           idx = ind3 * inmap.nx * inmap.ny + ind2 * inmap.nx + ind1;
           if (inmap.data[idx] >= threshold) {
               // Compute bead coordinate with a small Gaussian offset.
               asu_x = inmap.apix_x * ind1 + inmap.origin_x + inmap.apix_x * (gaussRand() - 0.5);
               asu_y = inmap.apix_y * ind2 + inmap.origin_y + inmap.apix_y * (gaussRand() - 0.5);
               asu_z = inmap.apix_z * ind3 + inmap.origin_z + inmap.apix_z * (gaussRand() - 0.5);
               // Check that the ASU bead is within the map boundaries.
               if (asu_x >= x_min && asu_x <= x_max &&
                   asu_y >= y_min && asu_y <= y_max &&
                   asu_z >= z_min && asu_z <= z_max)
                   valid = 1;
           }
       }
       // Place the ASU bead.
       xcoo[bead_count * 3 + 0] = asu_x;
       xcoo[bead_count * 3 + 1] = asu_y;
       xcoo[bead_count * 3 + 2] = asu_z;
       bead_count++;
   
       // If symmetry is enabled, try to add symmetric copies.
       if (use_symmetry) {
           int j = 1;
           while (j <= j_max && bead_count < natoms) {
               // Upward copy: rotate by +j*twist about (center_x, center_y) and shift z by +j*rise.
               double theta = j * twist_rad;
               double x_new = cos(theta) * (asu_x - center_x) - sin(theta) * (asu_y - center_y) + center_x;
               double y_new = sin(theta) * (asu_x - center_x) + cos(theta) * (asu_y - center_y) + center_y;
               double z_new = asu_z + j * rise;
               if (x_new >= x_min && x_new <= x_max &&
                   y_new >= y_min && y_new <= y_max &&
                   z_new >= z_min && z_new <= z_max) {
                   xcoo[bead_count * 3 + 0] = x_new;
                   xcoo[bead_count * 3 + 1] = y_new;
                   xcoo[bead_count * 3 + 2] = z_new;
                   bead_count++;
                   if (bead_count >= natoms)
                       break;
               }
               // Downward copy: rotate by -j*twist and shift z by -j*rise.
               theta = -j * twist_rad;
               x_new = cos(theta) * (asu_x - center_x) - sin(theta) * (asu_y - center_y) + center_x;
               y_new = sin(theta) * (asu_x - center_x) + cos(theta) * (asu_y - center_y) + center_y;
               z_new = asu_z - j * rise;
               if (x_new >= x_min && x_new <= x_max &&
                   y_new >= y_min && y_new <= y_max &&
                   z_new >= z_min && z_new <= z_max) {
                   xcoo[bead_count * 3 + 0] = x_new;
                   xcoo[bead_count * 3 + 1] = y_new;
                   xcoo[bead_count * 3 + 2] = z_new;
                   bead_count++;
               }
               j++;
           }
       }
   }

    // Write out the bead model in PDB format.
    //writePDB("beadmodel.pdb", xcoo, natoms);

    
    // Render bead model into a new density map.
    DensityMap beadMap;
    beadMap.nx = inmap.nx; beadMap.ny = inmap.ny; beadMap.nz = inmap.nz;
    beadMap.apix_x = inmap.apix_x; beadMap.apix_y = inmap.apix_y; beadMap.apix_z = inmap.apix_z;
    beadMap.origin_x = inmap.origin_x; beadMap.origin_y = inmap.origin_y; beadMap.origin_z = inmap.origin_z;
    beadMap.n = inmap.n;
    beadMap.data = (double *) calloc(beadMap.n, sizeof(double));
    if (!beadMap.data) { fprintf(stderr, "Memory allocation error\n"); exit(1); }
    //double sigma = inmap.apix_x / 2.0; // Set the Gaussian width.
    //renderBeadModel(&beadMap, natoms, xcoo, sigma);


    double *render_factor;
    // render_factor could be used in the future to scale bead density
    render_factor = (double *) calloc(natoms, sizeof(double));
    for (size_t i = 0; i < natoms; i++)  render_factor[i] = 1.0; 
    double peng_kernel_width = 10.0;  // adjust as needed (in Å)
    double peng_sigma = 2.0;          // adjust as needed
    renderPeng(&beadMap, peng_kernel_width, peng_sigma, natoms, atom_type, xcoo, render_factor);

        
    // Write out the bead map for debugging
    /*
    if (writeMRC("beadmap.mrc", &beadMap, header)) {
        fprintf(stderr, "Error writing output MRC file %s\n", outmap_fname);
        exit(1);
    } */

    // Match Structure Factor
    matchStructureFactor(&inmap, &beadMap, 100);

    // Match Density Histogram
    matchHistogram(&inmap, &beadMap);

    // Apply a smooth low-pass filter to the input map.
    smoothLowpassFilter_many(&inmap, resolution);

    // For this first version, we simply write out the processed input map.
    if (writeMRC(outmap_fname, &inmap, header)) {
        fprintf(stderr, "Error writing output MRC file %s\n", outmap_fname);
        exit(1);
    }
    
    // Clean up.
    free(inmap.data);
    free(beadMap.data);
    free(xcoo);
    
    return 0;
}

