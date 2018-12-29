![screenshot](cells.png)

# Celadro: Cells as active droplets

Phase-field modelling of epithelial cells using finite-difference integrator.

## Publications

TBD

## Building

We use cmake. Typically you would type:
```
mkdir build
cd build
cmake ..
make
```

We rely on the `boost::program_options` which must be installed prior to
building the program. We use modern C++17 features, such that you will
require a modern compiler.

## Input and running

The input of a simulation is a directory with the following structure:

```
./config.dat    <- configuration of the simulation
./cells/cell_type1.dat    <- configuration of the different cell types
./cells/cell_type2.dat
...
```

The names of the cell types files are arbitrary, they just need to have the
`.dat` extension and be located in the `cells/` subdirectory. See examples
that can be found in `examples/`.

The code is run from the command line as follows:

`./celadro input_dir [options]`

and produces output in `input_dir/data/`. Type `./celadro -h` for a list of
available options. The program supports compressed (using zip) output with the
option flag `--compression` or `-c`.

## Examples

Examples runs and ploting scripts can be found in the `examples/` directory. Try
to run `bash make_all_movies.sh` and have a look at the results.

## Multi-threading and CUDA

The code supports both multi-threading and CUDA.

The CUDA-enabled version of the code can be built using
```
mkdir build-cuda
cd build-cuda
cmake ..
make celadro-cuda
```
