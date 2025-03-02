CC = gcc
CFLAGS = -O2 -Wall -std=c99
LIBS = -lfftw3 -lm

TARGET = visdem

all: $(TARGET)

$(TARGET): visdem.c
	$(CC) $(CFLAGS) -o $(TARGET) visdem.c $(LIBS)

clean:
	rm -f $(TARGET)
