# © 2020 ETH Zurich, Mechanics and Materials Lab
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

add_executable(${PROJECT_NAME}-ElementsTests
               Bar_Test.cc
               CoreElement_Test.cc
               TwoNodeCorotationalBeam_Test.cc
               ForceElement_Test.cc
               TimoshenkoBeamElement_Test.cc
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
               mesh/refine_segment_mesh_Test.cc
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

find_package(GTest 1.8.1 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-ElementsTests PUBLIC
                      ${PROJECT_NAME}-elements
                      GTest::gmock_main
)
