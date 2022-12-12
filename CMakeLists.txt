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

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(ae108 VERSION 0.1.0 LANGUAGES C CXX)
enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

include(cmake/ae108_install_library.cmake)
include(cmake/ae108_add_cppptest.cmake)

add_subdirectory(cppptest)
foreach(AE108_LIBRARY elements cpppetsc cppslepc assembly solve cmdline)
    add_subdirectory(${AE108_LIBRARY})
endforeach()
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
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-program-options-dev,libeigen3-dev,librange-v3-dev,libxrender1,libvtk9-dev,libpetsc-${AE108_PETSC_SCALAR_TYPE}-dev,libslepc-${AE108_PETSC_SCALAR_TYPE}3.15-dev")
set(CPACK_PACKAGE_DESCRIPTION "C++ foundation for computational solid mechanics simulations using a variational framework, primarily focusing on the Finite Element Method (FEM)")
set(CPACK_PACKAGE_CONTACT "Mechanics & Materials Lab at ETH Zurich <gerhard.braeunlich@id.ethz.ch>")

include(CPack)
