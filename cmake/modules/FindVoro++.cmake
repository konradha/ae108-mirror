# © 2021 ETH Zurich, Mechanics and Materials Lab
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

if (NOT TARGET external::voro++)
  include(FindPackageHandleStandardArgs)
  if (NOT DEFINED ENV{Voro++_LIBRARIES} OR NOT DEFINED ENV{Voro++_INCLUDE_DIRS})
    find_library(Voro++_LIBRARIES NAMES voro++ PATHS $ENV{VOROPP_DIR}/lib)
    find_path(Voro++_INCLUDE_DIRS voro++.hh PATH_SUFFIXES voro++ PATHS $ENV{VOROPP_DIR}/include/voro++)
  else()
    find_library(Voro++_LIBRARIES NAMES voro++ PATHS $ENV{Voro++_LIBRARIES} NO_DEFAULT_PATH)
    find_path(Voro++_INCLUDE_DIRS voro++.hh PATH_SUFFIXES voro++ PATHS $ENV{Voro++_INCLUDE_DIRS} NO_DEFAULT_PATH)
  endif()

  find_package_handle_standard_args(Voro++ DEFAULT_MSG Voro++_LIBRARIES Voro++_INCLUDE_DIRS)

  if (Voro++_FOUND)
    add_library(external::voro++ UNKNOWN IMPORTED)
    set_property(TARGET external::voro++ PROPERTY IMPORTED_LOCATION ${Voro++_LIBRARIES})
    set_property(TARGET external::voro++ PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Voro++_INCLUDE_DIRS})
  endif()
endif()
