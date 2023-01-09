# © 2021 ETH Zurich, Mechanics and Materials Lab
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

add_executable(${PROJECT_NAME}-MeshingTests
        cppgmsh/construct_cylinders_Test.cc
        cppgmsh/construct_polyhedron_Test.cc
        cppgmsh/construct_rectangle_Test.cc
        cppgmsh/construct_box_Test.cc
        cppgmsh/copy_entities_Test.cc
        cppgmsh/extrude_entities_Test.cc
        cppgmsh/extrude_surface_Test.cc
        cppgmsh/fragment_entities_Test.cc
        cppgmsh/fuse_entities_Test.cc
        cppgmsh/intersect_entities_Test.cc
        cppgmsh/remove_entities_Test.cc
        cppgmsh/rotate_entities_Test.cc
        cppgmsh/translate_entities_Test.cc
        cppgmsh/generate_mesh_Test.cc
        cppgmsh/extract_mesh_Test.cc
        # cppgmsh/heal_periodic_surfaces_Test.cc
        # cppgmsh/heal_periodic_domains_Test.cc
        cppgmsh/as_affine_transform_Test.cc
        cppgmsh/get_boundary_of_Test.cc
        cppgmsh/get_centroid_of_Test.cc
        cppgmsh/get_coords_of_Test.cc
        cppgmsh/get_elements_in_Test.cc
        cppgmsh/get_entities_in_Test.cc
        cppgmsh/get_entities_of_Test.cc
        cppgmsh/get_nodes_of_Test.cc
        cppgmsh/get_normal_of_Test.cc
        cppgmsh/get_periodic_nodes_of_Test.cc
        cppgmsh/get_physical_groups_Test.cc
        cppgmsh/get_points_of_Test.cc
        cppgmsh/set_domain_entities_periodic_Test.cc
        cppgmsh/set_expert_mode_Test.cc
        cppgmsh/set_granularity_Test.cc
        cppgmsh/set_periodic_Test.cc
        cppgmsh/set_physical_group_of_Test.cc
        cppgmsh/synchronize_Test.cc
        cppgmsh/import_shapes_Test.cc
        cppgmsh/read_file_Test.cc
        cppgmsh/write_file_Test.cc
        BoundaryRepresentation_Test.cc
        BoundingBox_Test.cc
        bounding_box_of_Test.cc
        construct_rectilinear_grid_Test.cc
        construct_periodic_point_cloud_Test.cc
        construct_voronoi_cell_Test.cc
)

add_test(NAME ${PROJECT_NAME}-MeshingTests
         COMMAND $<TARGET_FILE:${PROJECT_NAME}-MeshingTests>)

find_package(GTest 1.8.1 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-MeshingTests
                      PRIVATE ${PROJECT_NAME}::meshing
                      PRIVATE GTest::gmock_main
)
