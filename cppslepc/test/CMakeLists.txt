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

add_executable(${PROJECT_NAME}-CppSlepcTests
        EigenvalueProblemSolverDivergedException_Test.cc
        computeEigenvalues_Test.cc
        main.cc
)

add_test(NAME ${PROJECT_NAME}-CppSlepcTests
         COMMAND $<TARGET_FILE:${PROJECT_NAME}-CppSlepcTests>)

find_package(AE108_GTest MODULE REQUIRED)
target_link_libraries(${PROJECT_NAME}-CppSlepcTests
                      PRIVATE ${PROJECT_NAME}::cppslepc
                      PRIVATE GTest::gmock
)