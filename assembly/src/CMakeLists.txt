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

add_library(${PROJECT_NAME}-assembly
            ./AssembleEnergyPlugin.cc
            ./AssembleForceIfPlugin.cc
            ./AssembleForceVectorPlugin.cc
            ./AssembleLumpedMassMatrixPlugin.cc
            ./AssembleConsistentMassMatrixPlugin.cc
            ./AssembleStiffnessMatrixPlugin.cc
            ./Assembler.cc
            ./AssemblerGroup.cc
            ./AssemblerTypeTraits.cc
            ./ConcatenateFeaturePlugins.cc
            ./DefaultFeaturePlugins.cc
            ./DerivePluginsUniquely.cc
            ./DeriveUniquely.cc
            ./FeaturePlugin.cc
            ./FeaturePlugins.cc
            ./HasUniqueTypeTrait.cc
            ./IsSingleUniqueType.cc
            ./StaticLooper.cc
            ./UpdateInternalVariablesPlugin.cc
            ./deserialize.cc
            ./groupElementDataPerVertex.cc
            ./resizeIfPossible.cc
            ./serialize.cc
)

add_library(${PROJECT_NAME}::assembly ALIAS ${PROJECT_NAME}-assembly)

target_include_directories(${PROJECT_NAME}-assembly PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

find_package(Eigen3 3.3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-assembly
                      PUBLIC ${PROJECT_NAME}::cpppetsc
                      PUBLIC Eigen3::Eigen
)

target_compile_features(${PROJECT_NAME}-assembly PUBLIC cxx_std_17)

ae108_install_library(assembly)
