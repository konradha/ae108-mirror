# © 2020 ETH Zurich, Mechanics and Materials Lab
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

project(ae108-examples LANGUAGES C CXX)

if (NOT TARGET ae108::elements)
    find_package(ae108 REQUIRED)
endif()

add_executable(${PROJECT_NAME}-Cmdline Cmdline.cc)
target_link_libraries(${PROJECT_NAME}-Cmdline
                      PRIVATE ae108::cmdline
)

add_executable(${PROJECT_NAME}-MeshGeneration MeshGeneration.cc)
target_link_libraries(${PROJECT_NAME}-MeshGeneration
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-Basic Basic.cc)
target_link_libraries(${PROJECT_NAME}-Basic
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-Output Output.cc)
target_link_libraries(${PROJECT_NAME}-Output
                      PRIVATE ae108::cpppetsc
)

add_executable(${PROJECT_NAME}-Input Input.cc)
target_link_libraries(${PROJECT_NAME}-Input
                      PRIVATE ae108::cpppetsc
)