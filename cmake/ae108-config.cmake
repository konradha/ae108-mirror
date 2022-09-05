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

include(CMakeFindDependencyMacro)

find_dependency(Boost 1.67 COMPONENTS program_options REQUIRED)

find_dependency(Eigen3 3.4 CONFIG REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
find_dependency(AE108_PETSc MODULE 3.15)

find_dependency(range-v3 0.11.0 CONFIG REQUIRED)

foreach(AE108_LIBRARY elements cpppetsc cppslepc assembly solve cmdline)
    include("${CMAKE_CURRENT_LIST_DIR}/ae108-${AE108_LIBRARY}-export.cmake")
endforeach()
