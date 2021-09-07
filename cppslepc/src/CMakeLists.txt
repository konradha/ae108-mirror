# © 2021 ETH Zurich, Mechanics and Materials Lab
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

add_library(${PROJECT_NAME}-cppslepc
    Context.cc
    LinearEigenvalueProblemSolver.cc
    computeEigenvalues.cc
)

target_include_directories(${PROJECT_NAME}-cppslepc PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(${PROJECT_NAME}-cppslepc
                        PUBLIC cxx_std_11
)

find_package(AE108_SLEPc 1.12 MODULE REQUIRED)
target_link_libraries(${PROJECT_NAME}-cppslepc
                      PUBLIC ${PROJECT_NAME}::external::slepc
                      PUBLIC ${PROJECT_NAME}::cpppetsc
)

add_library(${PROJECT_NAME}::cppslepc ALIAS ${PROJECT_NAME}-cppslepc)

ae108_install_library(cppslepc)
