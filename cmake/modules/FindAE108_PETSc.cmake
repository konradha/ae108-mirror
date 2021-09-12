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

    macro(try_compile_with_petsc COMPILE_DEFINITIONS)
        try_compile(
            AE108_PETSc_COMPILE_RESULT
            "${CMAKE_CURRENT_BINARY_DIR}/petsc"
            SOURCES "${CMAKE_CURRENT_LIST_DIR}/AE108_PETSc.cc"
            COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS}"
            CXX_STANDARD 11
            CXX_STANDARD_REQUIRED TRUE
            LINK_LIBRARIES ae108::external::petsc
        )
    endmacro()
    
    try_compile_with_petsc("")
    if(NOT AE108_PETSc_COMPILE_RESULT)
        message(FATAL_ERROR "Compiling a simple executable that uses PETSc failed.")
    endif()

    try_compile_with_petsc("-DAE108_PETSC_COMPLEX")
    if(AE108_PETSc_COMPILE_RESULT)
        target_compile_definitions(ae108::external::petsc
            INTERFACE AE108_PETSC_COMPLEX=1
        )
        message(STATUS "A version of PETSc with complex scalar type was found.")
    else()
        try_compile_with_petsc("-DAE108_PETSC_REAL")
        if(NOT AE108_PETSc_COMPILE_RESULT)
            MESSAGE(WARNING "PETSc uses an incompatible scalar type.")
        endif()
    endif()
endif()
