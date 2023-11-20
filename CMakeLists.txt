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

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(ae108 VERSION 0.1.0 LANGUAGES C CXX)
enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

include(cmake/ae108_install_library.cmake)
include(cmake/ae108_add_cppptest.cmake)

option(USE_MESHING "Build meshing module" ON)

add_subdirectory(cppptest)
foreach(AE108_LIBRARY elements cpppetsc cppslepc assembly solve cmdline)
    add_subdirectory(${AE108_LIBRARY})
endforeach()
if (${USE_MESHING})
    add_subdirectory(meshing)
endif()
add_subdirectory(examples)

add_test(NAME ${PROJECT_NAME}-ExamplesTests COMMAND "${CMAKE_CURRENT_LIST_DIR}/tests/run.py")

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
write_basic_package_version_file(${PROJECT_NAME}-config-version.cmake
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion
)
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
        cmake/${PROJECT_NAME}-config.cmake
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

install(FILES
        cmake/modules/FindAE108_PETSc.cmake
        cmake/modules/FindAE108_SLEPc.cmake
        cmake/modules/AE108_PETSc.cc
        cmake/modules/FindVoro++.cmake
        cmake/modules/FindGmsh.cmake
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/modules"
)

install(FILES LICENSE DESTINATION .)

if(${AE108_PETSC_COMPLEX})
  set(AE108_PETSC_SCALAR_TYPE "complex")
else()
  set(AE108_PETSC_SCALAR_TYPE "real")
endif()

set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_NAME "libae108-${AE108_PETSC_SCALAR_TYPE}-dev")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-program-options-dev,libeigen3-dev,librange-v3-dev,libxrender1,libvtk9-dev,libpetsc-${AE108_PETSC_SCALAR_TYPE}-dev,libslepc-${AE108_PETSC_SCALAR_TYPE}3.15-dev,libgmsh-dev,voro++-dev")
set(CPACK_PACKAGE_DESCRIPTION "C++ foundation for computational solid mechanics simulations using a variational framework, primarily focusing on the Finite Element Method (FEM)")
set(CPACK_PACKAGE_CONTACT "Mechanics & Materials Lab at ETH Zurich <gerhard.braeunlich@id.ethz.ch>")

include(CPack)
