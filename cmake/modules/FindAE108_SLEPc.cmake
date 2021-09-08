# © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
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
pkg_check_modules(SLEPc SLEPc IMPORTED_TARGET)
find_package(AE108_PETSc 1.12 MODULE)

find_package_handle_standard_args(AE108_SLEPc
                                  REQUIRED_VARS SLEPc_FOUND AE108_PETSc_FOUND
                                  VERSION_VAR SLEPc_VERSION
)

if(AE108_SLEPc_FOUND AND NOT TARGET ae108::external::slepc)
    add_library(ae108::external::slepc INTERFACE IMPORTED)
    target_link_libraries(ae108::external::slepc
                          INTERFACE PkgConfig::SLEPc
                          INTERFACE ae108::external::petsc
    )
endif()