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

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

add_library(${PROJECT_NAME}-cppslepc
    Context.cc
    EigenPair.cc
    EigenvalueSolverDivergedException.cc
    InvalidEigenvalueIndexException.cc
    InvalidProblemTypeException.cc
    LinearEigenvalueProblemSolver.cc
    NoOperatorsSetException.cc
    computeEigenvalues.cc
    computeSmallestEigenvalues.cc
)

target_include_directories(${PROJECT_NAME}-cppslepc PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(${PROJECT_NAME}-cppslepc
                        PUBLIC cxx_std_17
)

find_package(AE108_SLEPc 3.15 MODULE REQUIRED)
target_link_libraries(${PROJECT_NAME}-cppslepc
                      PUBLIC ${PROJECT_NAME}::external::slepc
                      PUBLIC ${PROJECT_NAME}::cpppetsc
)

add_library(${PROJECT_NAME}::cppslepc ALIAS ${PROJECT_NAME}-cppslepc)

ae108_install_library(cppslepc)
