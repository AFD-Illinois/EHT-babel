# EHT-babel
Tools for converting output from ```KORAL```, ```bhlight```, and ```BHAC``` to the HARM HDF5 format specified [here](https://github.com/AFD-Illinois/docs/wiki/Output-Format).

## Usage

### BHAC:

Currently supports only BHAC binary-format files using the ```LOG_KS``` coordinate system.

Usage:

```bash
$ python3 bhac-translate.py [--gridfile] file1.blk [file2.blk ...]
```

Outputs HARM ```.h5``` files alongside all BHAC ```.blk``` files specified.  If passed ```--gridfile```, also outputs a file ```grid.h5``` containing the locations of grid zone centers in KS coordinates, and the local metric in native coordinates.

### KORAL:

```bash
./babel path/to/define.h path/to/koral.dat
```

You may need to comment/uncomment the type of KORAL data you expect to read. Do this by modifying the lines in ```babel.c```

```c
geom.dtype = KORAL_GRMHD;
geom.dtype = KORAL_RADGRMHD;
```

### bhlight2d:

Currently supports only ASCII-format dump files using MKS coordinates.

Usage:

```bash
$ python3 bhlight2d-translate.py file1 [file2 ...]
```

Outputs HARM ```.h5``` files alongside each specified text file.
