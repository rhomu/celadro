![screenshot](cells.png)

# Celadro: Cells as active droplets

Phase-field modelling of epithelial cells using finite-difference integrator.<br/>
[Example movie 1](assets/movie1.mp4?raw=true) [Example movie 2](assets/movie2.mp4?raw=true)

## Publications

This code has served as a basis for the simulations included in the following papers:

- Active Inter-cellular Forces in Collective Cell Motility  
  G Zhang, R Mueller, A Doostmohammadi, JM Yeomans  
  arXiv preprint arXiv:2005.13087, 2020
- Emergence of active nematic behavior in monolayers of isotropic cells  
  R Mueller, JM Yeomans, A Doostmohammadi  
  Physical review letters 122 (4), 048004  
  [Movie 1](https://arxiv.org/src/1811.05040v2/anc/movie_1.mp4) [Movie 1](https://arxiv.org/src/1811.05040v2/anc/movie_2.mp4)
- Sustained oscillations of epithelial cell sheets  
  G Peyret, R Mueller, J d'Alessandro, et al.  
  Biophysical journal 117 (3), 464-478

## Building

We use cmake. Typically you would type:
```
mkdir build
cd build
cmake ..
make
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
flag `--compression` or `-c`. Typical usage is `./celadro runcard.dat -fco output/`

Type `./celadro -h` for a list of available options.

## Examples

Examples runs and ploting scripts can be found in the `example` directory. Try
to run `bash make_all_movies.sh` and have a look at the results.

## Multi-threading and CUDA

WARNING !!! this feature is currently unmaintained and can not be used for now (2020) !!!

The code supports both multi-threading and CUDA. The CUDA-enabled version of
the code can be built using
```
mkdir build-cuda
cd build-cuda
cmake ..
make celadro-cuda
```

## Troubleshooting

### Hard to solve linking problems with `boost::program_options`

On some platforms the compiler might be unable to link to the boost library properly
which results in `undefined reference` errors. This can be due to the installed boost
libraries using a different ABI than the code when it is compiled with a recent
compiler. In order to solve this problem, either switch back to an older compiler
(e.g. `gcc-4.9`) or uncomment the line
```
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
```
in `CMakeLists.txt`.
