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

add_executable(${PROJECT_NAME}-AssemblyTests
               ./AssemblerGroup_Test.cc
               ./Assembler_Test.cc
               ./ConcatenateFeaturePlugins_Test.cc
               ./DerivePluginsUniquely_Test.cc
               ./DeriveUniquely_Test.cc
               ./FeaturePlugin_Test.cc
               ./HasUniqueTypeTrait_Test.cc
               ./IsSingleTypeTrait_Test.cc
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
