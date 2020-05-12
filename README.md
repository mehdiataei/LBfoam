![Image of a 2D foam structure made by LBfoam](imgs/2Dfoam.png)
# LBfoam

LBfoam is an open-source CFD solver based on the lattice Boltzmann method for foaming simulations. The solver is an extended version of the [Palabos](https://palabos.unige.ch/) library.

# Installation

LBfoam installation is very similar to the Palabos library and it does not depend on any external dependencies. LBfoam uses [`Scons`](https://scons.org/) build tool.

The mandatory packages for installation are `gcc` (or `clang`), `make`, `python3`. For MPI parallel computations, `libopenmpi` library is required. To output results in `.gif` format, `imagemagick` library must be installed.

For Debian based distributions, the following command can install the required libraries.

```
$ sudo apt install gcc python3 make imagemagick libopenmpi-dev
```


Similar to Palabos the output of the simulations are in the form of `VTK` files, which can be visualized using the
of [Paraview](https://www.paraview.org/) software.

## Run an example

First, close the repository using the following git command:

```
$ git clone https://github.com/mehdiataei/LBfoam.git
```


Compile the `bucket2D` example:
```
$ cd examples/lbfoam/bucket2D
$ make
```

(Note: To compile the software on MacOS, uncomment the ` -DPLB_MAC_OS_X` compilation flag in the Makefile).

Run the example using the following command:

```
$ ./bucket2D bucket2D.xml

```

To run the example in parallel using 8 cores for example:

```
$ ./mpirun -np 8 bucket2D bucket2D.xml

```

# Getting help and bug report

Please submit an issue if you found a bug in the program or needed help with the software.