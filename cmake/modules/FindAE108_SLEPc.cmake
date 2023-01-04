# © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
# © 2020 California Institute of Technology
#
# This file is part of ae108.
#
# ae108 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any
# later version.
#
# ae108 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ae108. If not, see <https://www.gnu.org/licenses/>.

include(FindPkgConfig)
pkg_check_modules(SLEPc slepc IMPORTED_TARGET)
find_package(AE108_PETSc 3.15 MODULE)

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
