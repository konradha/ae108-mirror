# © 2020 ETH Zurich, Mechanics and Materials Lab
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

include(CMakeFindDependencyMacro)

find_dependency(Boost 1.67 COMPONENTS program_options REQUIRED)

find_dependency(Eigen3 3.4 CONFIG REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
find_dependency(AE108_PETSc MODULE 3.15)
find_dependency(AE108_SLEPc MODULE 3.15)

find_dependency(range-v3 0.11.0 CONFIG REQUIRED)

foreach(AE108_LIBRARY elements cpppetsc cppslepc assembly solve cmdline)
    include("${CMAKE_CURRENT_LIST_DIR}/ae108-${AE108_LIBRARY}-export.cmake")
endforeach()
