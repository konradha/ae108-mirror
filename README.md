<!---
 © 2020 ETH Zurich, Mechanics and Materials Lab
 © 2020 California Institute of Technology

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
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

- ```ae108-assembly```: Assemble local data.

    A library for assembling MPI-local data in ```cpppetsc``` meshes. For instance, assemble the local energy by summing the energy of all local element instances.

- ```ae108-solve```: Minimize the total energy.

    A library for minimizing the assembled total energy using PETSc's solvers or optimizers as provided by ```cpppetsc```.

## Installation

The project uses [CMake](https://cmake.org) as its build system generator. The following third party libraries are required and located using CMake's ```find_package```.

- [Boost](https://www.boost.org) (component program_options): version 1.67
- [Eigen](http://eigen.tuxfamily.org): version 3.3
- [Google Test](https://github.com/google/googletest): version 1.8.1
- [MPI](https://cmake.org/cmake/help/latest/module/FindMPI.html): version 3.1
- [PETSc](https://www.mcs.anl.gov/petsc/): version 3.10

Of course, these libraries are covered by their own license terms. Since PETSc does not provide a CMake configuration file, PETSc is found using the provided find module in ```cmake/modules/```, which in turn is based on ```pkg-config```.

Once you have installed these libraries, run CMake to build the project, choosing a location to install the library to by specifying ```CMAKE_INSTALL_PREFIX```; see the following example. Of course, depending on your setup, you might need to add a ```-DCMAKE_PREFIX_PATH='...'``` parameter to tell CMake the location of the third party library installations.

```bash
mkdir build
cd build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX='install/to/path' \
    ..
```

Now run ```make``` to build and install:

```bash
make -j$(nproc) install
```

## Usage

To link and use the libraries point CMake to the installation using ```-DCMAKE_PREFIX_PATH='...'``` and use CMake's ```find_package``` to find the package ```ae108```. Now link the libraries that you want to link using CMake's ```target_link_libraries```. The following targets are installed:

- ```ae108::cmdline```
- ```ae108::cpppetsc```
- ```ae108::assembly```
- ```ae108::solve```

## License

AE108 is available under the [Apache 2.0 license](https://choosealicense.com/licenses/apache-2.0/).

## Releases

The latest release is v0.1.0. It is also available via the ETH library: [ae108](https://search.library.ethz.ch/permalink/f/13kse66/data_archiveIE15605648).
