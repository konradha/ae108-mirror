# © 2020 ETH Zurich, Mechanics and Materials Lab
# © 2020 California Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

include(FindPkgConfig)
pkg_check_modules(PETSc PETSc IMPORTED_TARGET)
find_package(MPI 3.1)
find_package(HDF5 1.10 MODULE)

find_package_handle_standard_args(AE108_PETSc
                                  REQUIRED_VARS PETSc_FOUND MPI_CXX_FOUND HDF5_FOUND
                                  VERSION_VAR PETSc_VERSION
)

if(AE108_PETSc_FOUND AND NOT TARGET ae108::external::petsc)
    add_library(ae108::external::petsc INTERFACE IMPORTED)
    target_include_directories(ae108::external::petsc
                               INTERFACE ${HDF5_INCLUDE_DIRS}
    )
    target_link_libraries(ae108::external::petsc
                          INTERFACE MPI::MPI_CXX
                          INTERFACE PkgConfig::PETSc
    )
endif()
