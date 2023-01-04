<!---
 © 2020 ETH Zurich, Mechanics and Materials Lab
 © 2020 California Institute of Technology

 This file is part of ae108.

 ae108 is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or any
 later version.

 ae108 is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ae108. If not, see <https://www.gnu.org/licenses/>.
-->

# AE108

## Introduction

The aim of this project is to provide a C++ foundation for computational solid mechanics simulations using a variational framework, primarily focusing on the Finite Element Method (FEM).
See the [Computational Solid Mechanics lecture notes](https://www.mm.ethz.ch/) for more information about the approach that is used.
The code is developed by the Mechanics & Materials Lab at ETH Zurich.

## Summary

The following libraries are provided:

- ```ae108-cmdline```: Read command line parameters.

    A simple library for reading command line parameters with few lines of code (based on Boost.Program_options).

- ```ae108-cpppetsc```: Create a mesh.

    A library for creating a FEM mesh and managing the degrees of freedom associated with the elements and their vertices (based on PETSc's DMPLEX). It also provides wrappers around PETSc's solvers (KSP, SNES, and TAO).

- ```ae108-elements```: Define elements.

    A library for specifying the behaviour of the elements that make up the mesh. For instance, use this library to define the energy needed to deform an element.

- ```ae108-assembly```: Assemble local data.

    A library for assembling MPI-local data in ```cpppetsc``` meshes. For instance, assemble the local energy by summing the energy of all local element instances.

- ```ae108-solve```: Minimize the total energy.

    A library for minimizing the assembled total energy using PETSc's solvers or optimizers as provided by ```cpppetsc```.
    
## Installation

🐧: [deb packages](/../../packages/) available for Ubuntu jammy.
The packages can be downloaded and installed by
`apt-get -f install ./package.deb`.

🪟 / 🍏: You have to compile it by yourself.

The project uses [CMake](https://cmake.org) as its build system generator. The following third party libraries are required and located using CMake's ```find_package```.

- [Boost](https://www.boost.org) (component program_options): version 1.67
- [Eigen](http://eigen.tuxfamily.org): version 3.4
- [Google Test](https://github.com/google/googletest): version 1.8.1
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/): version 1.10
- [MPI](https://cmake.org/cmake/help/latest/module/FindMPI.html): version 3.1
- [PETSc](https://www.mcs.anl.gov/petsc/): version 3.15
- [Range-v3](https://github.com/ericniebler/range-v3): version 0.11
- [SLEPc](https://slepc.upv.es/): version 3.15

Of course, these libraries are covered by their own license terms. Since PETSc and SLEPc do not provide a CMake configuration file, these libraries are found using the provided find modules in ```cmake/modules/```, which in turn are based on ```pkg-config```.

Once you have installed these libraries, run CMake to build the project, choosing a location to install the library to by specifying ```CMAKE_INSTALL_PREFIX```; see the following example. Of course, depending on your setup, you might need to add a ```-DCMAKE_PREFIX_PATH='...'``` parameter to tell CMake the location of the third party library installations.

```bash
cmake \
    -S . \
    -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX='install/to/path'
```

Now use ```cmake``` to build and install:

```bash
cmake --build build --target install -- -j$(nproc)
```

## Usage

To link and use the libraries point CMake to the installation using ```-DCMAKE_PREFIX_PATH='...'``` and use CMake's ```find_package``` to find the package ```ae108```. Now link the libraries that you want to link using CMake's ```target_link_libraries```. The following targets are installed:

- ```ae108::cmdline```
- ```ae108::cpppetsc```
- ```ae108::elements```
- ```ae108::assembly```
- ```ae108::solve```

## License

AE108 is available under the [GNU General Public License v3.0 or later](https://www.gnu.org/licenses/gpl-3.0.html).

## Releases

The releases of this project are are available [here](https://gitlab.ethz.ch/mechanics-and-materials/ae108/-/releases). Release v0.1.0 is also available via the ETH library: [ae108](https://search.library.ethz.ch/permalink/f/13kse66/data_archiveIE15605648).
