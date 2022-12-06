.. raw:: html

   <!---
    © 2022 ETH Zurich, Mechanics and Materials Lab

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

Libraries
=========

The following libraries are provided as part of *ae108*:

-  ``ae108-cmdline``: Read command line parameters.

   A simple library for reading command line parameters with few lines
   of code (based on Boost.Program_options).
   
   
-  ``ae108-cpppetsc``: Manage a distributed data mesh.

   A library for creating a FEM mesh and managing the degrees of freedom
   associated with the elements and their vertices (based on PETSc’s
   DMPLEX). It also provides wrappers around PETSc’s solvers (KSP, SNES,
   and TAO).

-  ``ae108-elements``: Define elements.

   A library for specifying the behaviour of the elements that make up
   the mesh. For instance, use this library to define the energy needed
   to deform an element.

-  ``ae108-assembly``: Assemble local data.

   A library for assembling MPI-local data in ``cpppetsc`` meshes. For
   instance, assemble the local energy by summing the energy of all
   local element instances.

-  ``ae108-solve``: Minimize the total energy.

   A library for minimizing the assembled total energy using PETSc’s
   solvers or optimizers as provided by ``cpppetsc``.


.. toctree::
   :maxdepth: 1

   ae108-cmdline <../cmdline/README.rst>
   ae108-elements <../elements/README.rst>