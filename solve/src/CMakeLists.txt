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

add_library(${PROJECT_NAME}-solve
            ./AffineTransform.cc
            ./DampingParameters.cc
            ./DynamicSolver.cc
            ./DynamicState.cc
            ./InconsistentBoundaryConditionsException.cc
            ./NewmarkParameters.cc
            ./NonlinearSolver.cc
            ./TAOSolver.cc
            ./TransformingSolver.cc
            ./boundaryConditionsToTransform.cc
)

add_library(${PROJECT_NAME}::solve ALIAS ${PROJECT_NAME}-solve)

target_include_directories(${PROJECT_NAME}-solve PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(${PROJECT_NAME}-solve PUBLIC ${PROJECT_NAME}::assembly)

target_compile_features(${PROJECT_NAME}-solve PUBLIC cxx_std_11)

ae108_install_library(solve)
