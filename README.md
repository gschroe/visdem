# VISDEM 1.2 (2025)

**VISDEM 1.2** is a cryo‑EM density map sharpening tool developed by Gunnar F. Schroeder, Michaela Spiegel, and Amudha Duraisamy. It estimates the underlying atomic distribution of a protein from the density map and computes the corresponding radially averaged structure factor and density histogram. These calculated profiles are then used to adjust the input map’s spectrum and histogram, enhancing map clarity. VISDEM 1.2 also supports the application of helical symmetry during the atomic distribution estimation.

When using VISDEM please cite: \
M. Spiegel, A.K. Duraisamy, and G.F. Schröder. Improving the Visualisation of Cryo-EM Density Reconstructions \
J.Struct.Biol. (2015), 191(2):207-213


## Requirements

- A C compiler (GCC is recommended)
- [FFTW3](http://www.fftw.org/) library
- GNU Make

## Building

Clone the repository and build the program by running:

```bash
make
```

This will produce an executable named `visdem`.

## Usage

```
./visdem --input <input_map.mrc> --mass <mass in kDa> --resolution <cutoff in Å> --output <output_map.mrc> [--threshold <value>] [--rise <value> --twist <value>]
```

### Command-Line Options

- `-i, --input`  
  Input MRC map file.

- `-o, --output`  
  Output MRC map file.

- `-m, --mass`  
  Particle mass in kilodaltons (kDa).

- `-r, --resolution`  
  Resolution cutoff in Å.

- `-t, --threshold`  
  (Optional) Density threshold.

- `--rise`  
  (Optional) Rise parameter for helical symmetry.

- `--twist`  
  (Optional) Twist parameter for helical symmetry.

### Example

To run the program with helical symmetry using a rise of 4.75 Å and a twist of -1.2°:

```bash
./visdem -i map.mrc -m 400 -r 2.5 -o v.mrc --rise 4.75 --twist -1.2
```

If no threshold is provided, the program will compute one automatically from volume based on the provided mass.

## License

This project is licensed under the GNU General Public License (GPL) version 3. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

VISDEM 1.2 (2025) is developed by Gunnar F. Schroeder, Michaela Spiegel, and Amudha Duraisamy.

