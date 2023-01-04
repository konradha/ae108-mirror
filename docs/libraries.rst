.. raw:: html

   <!---
    © 2022 ETH Zurich, Mechanics and Materials Lab

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