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

add_executable(${PROJECT_NAME}-ElementsTests
               CoreElement_Test.cc
               ForceElement_Test.cc
               TimoshenkoBeamElement_Test.cc
               TimoshenkoBeamElementWithMass_Test.cc
               Minimal_Test.cc
               embedding/IsoparametricEmbedding_Test.cc
               integrator/IsoperimetricIntegrator_Test.cc
               materialmodels/Hookean_Test.cc
               materialmodels/Minimal_Test.cc
               mesh/Mesh_Test.cc
               mesh/generate_cuboid_mesh_Test.cc
               mesh/generate_quadratic_tetrahedron_mesh_Test.cc
               mesh/generate_quadratic_triangle_mesh_Test.cc
               mesh/generate_tetrahedron_mesh_Test.cc
               mesh/generate_triangle_mesh_Test.cc
               quadrature/Quadrature_Test.cc
               shape/Hexa8_Test.cc
               shape/Quad4_Test.cc
               shape/Seg2_Test.cc
               shape/Tet10_Test.cc
               shape/Tet4_Test.cc
               shape/Tri3_Test.cc
               shape/Tri6_Test.cc
               tensor/Tensor_Test.cc
               tensor/as_matrix_of_columns_Test.cc
               tensor/as_matrix_of_rows_Test.cc
               tensor/as_two_tensor_Test.cc
               tensor/as_vector_Test.cc
)

add_test(NAME ${PROJECT_NAME}-ElementsTests COMMAND $<TARGET_FILE:${PROJECT_NAME}-ElementsTests>)

find_package(AE108_GTest MODULE REQUIRED)
target_link_libraries(${PROJECT_NAME}-ElementsTests PUBLIC
                      ${PROJECT_NAME}-elements
                      GTest::gmock_main
)
