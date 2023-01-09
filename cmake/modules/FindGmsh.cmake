# © 2021 ETH Zurich, Mechanics and Materials Lab
#
# Licensed under the GNU General Public License (GPL), Version 2.0 (the
# "License"); you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# https:#www.gnu.org/licenses/old-licenses/gpl-2.0.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

if (NOT TARGET external::gmsh)
  include(FindPackageHandleStandardArgs)

  if (NOT DEFINED ENV{GMSH_LIBRARIES} OR NOT DEFINED ENV{GMSH_INCLUDE_DIRS})
    find_library(Gmsh_LIBRARIES NAMES gmsh)
    find_path(Gmsh_INCLUDE_DIRS gmsh.h)
  else()
    find_library(Gmsh_LIBRARIES NAMES gmsh PATHS $ENV{GMSH_LIBRARIES} NO_DEFAULT_PATH)
    find_path(Gmsh_INCLUDE_DIRS gmsh.h PATHS $ENV{GMSH_INCLUDE_DIRS} NO_DEFAULT_PATH)
  endif()

  find_package_handle_standard_args(Gmsh Gmsh_LIBRARIES Gmsh_INCLUDE_DIRS)

  if (Gmsh_FOUND)
    add_library(external::gmsh UNKNOWN IMPORTED)
    set_property(TARGET external::gmsh PROPERTY IMPORTED_LOCATION ${Gmsh_LIBRARIES})
    set_property(TARGET external::gmsh PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Gmsh_INCLUDE_DIRS})
  endif()
endif()


