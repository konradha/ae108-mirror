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

add_library(${PROJECT_NAME}-cppptest STATIC
            ./Matchers.cc
            ./isLocal.cc
            ./main.cc
)

add_library(${PROJECT_NAME}::cppptest ALIAS ${PROJECT_NAME}-cppptest)

target_include_directories(${PROJECT_NAME}-cppptest PUBLIC include)

find_package(AE108_PETSc 3.10 MODULE REQUIRED)
find_package(GTest 1.8.1 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-cppptest
                      PRIVATE ${PROJECT_NAME}::external::petsc
                      PUBLIC GTest::gmock
)

target_compile_features(${PROJECT_NAME}-cppptest PUBLIC cxx_std_11)
