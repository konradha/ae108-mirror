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

add_executable(${PROJECT_NAME}-AssemblyTests
               ./AssemblerGroup_Test.cc
               ./Assembler_Test.cc
               ./ConcatenateFeaturePlugins_Test.cc
               ./DerivePluginsUniquely_Test.cc
               ./DeriveUniquely_Test.cc
               ./HasUniqueTypeTrait_Test.cc
               ./IsSingleTypeTrait_Test.cc
               ./StaticLooper_Test.cc
               ./deserialize_Test.cc
               ./resizeIfPossible_Test.cc
               ./serialize_Test.cc
)

target_include_directories(${PROJECT_NAME}-AssemblyTests PRIVATE include)

target_link_libraries(${PROJECT_NAME}-AssemblyTests
                      PRIVATE ${PROJECT_NAME}::assembly
                      PRIVATE ${PROJECT_NAME}::cppptest
)

a108_add_cppptest(${PROJECT_NAME}-AssemblyTests)
