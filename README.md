# Celadro: Cells as active droplets

Phase-field modelling of epithelial cells using finite-difference integrator.

## Description

TBD in an upcoming publication.

## Building

We use cmake. Typically you would type:
```
mkdir build
cd build
cmake ..
```

We rely on the `boost::program_options` which must be installed prior to
building the program. We also use modern C++ features, such that you will
require a modern compiler (tested with g++-4.9).

## Running

The code is run from the command line and a runcard must always be given as the
first argument:

`./celadro runcard.dat [options]`

A runcard is a simple file providing the parameters for the run. Example
runcards can be found in the `example/` directory. Every option can also be
given to the program using the command line as `./celadro runcard.dat --option=arg`.
A complete list of available options can be obtained by typing `./celadro -h`.

By default the program writes output files in the current directory. This can be
changed using `--output=dir/` or `-o dir/`, where `dir/` is the target
directory. The program also supports compressed (using zip) output with the option
flag `--compression` or `-c`. When compression is on, the output name does not
denote the target directory but rather the file name of the compressed archive.
Note that with the current implementation, compression is not recommended for
long runs as the time to add a file to the archive grows with the size of the
archive.

Type `./celadro -h` for a list of available options.
