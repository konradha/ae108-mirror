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

add_library(${PROJECT_NAME}-assembly
            ./AssembleEnergyPlugin.cc
            ./AssembleForceIfPlugin.cc
            ./AssembleForceVectorPlugin.cc
            ./AssembleLumpedMassMatrixPlugin.cc
            ./AssembleConsistentMassMatrixPlugin.cc
            ./AssembleMassMatrixPlugin.cc
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

find_package(Eigen3 3.4 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-assembly
                      PUBLIC ${PROJECT_NAME}::cpppetsc
                      PUBLIC Eigen3::Eigen
)

target_compile_features(${PROJECT_NAME}-assembly PUBLIC cxx_std_17)

ae108_install_library(assembly)
