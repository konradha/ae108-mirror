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

add_library(${PROJECT_NAME}-meshing
    cppgmsh/Context.cc
    cppgmsh/Node.cc
    cppgmsh/as_affine_transform.cc
    cppgmsh/construct_box.cc
    cppgmsh/construct_cylinders.cc
    cppgmsh/construct_polyhedron.cc
    cppgmsh/construct_rectangle.cc
    cppgmsh/copy_entities.cc
    cppgmsh/extrude_entities.cc
    cppgmsh/extrude_surface.cc
    cppgmsh/fuse_entities.cc
    cppgmsh/fragment_entities.cc
    cppgmsh/intersect_entities.cc
    cppgmsh/remove_entities.cc
    cppgmsh/rotate_entities.cc
    cppgmsh/translate_entities.cc
    cppgmsh/get_boundary_of.cc
    cppgmsh/get_centroid_of.cc
    cppgmsh/get_coords_of.cc
    cppgmsh/get_elements_in.cc
    cppgmsh/get_entities_in.cc
    cppgmsh/get_entities_of.cc
    cppgmsh/get_nodes_of.cc
    cppgmsh/get_normal_of.cc
    cppgmsh/get_periodic_nodes_of.cc
    cppgmsh/get_physical_groups.cc
    cppgmsh/get_points_of.cc
    cppgmsh/set_domain_entities_periodic.cc
    cppgmsh/set_expert_mode.cc
    cppgmsh/set_granularity.cc
    cppgmsh/set_periodic.cc
    cppgmsh/set_physical_group_of.cc
    cppgmsh/heal_periodic_surfaces.cc
    cppgmsh/heal_periodic_domains.cc
    cppgmsh/import_shapes.cc
    cppgmsh/read_file.cc
    cppgmsh/write_file.cc
    cppgmsh/extract_mesh.cc
    cppgmsh/generate_mesh.cc
    cppgmsh/synchronize.cc
    BoundaryRepresentation.cc
    BoundingBox.cc
    bounding_box_of.cc
    construct_rectilinear_grid.cc
    construct_periodic_point_cloud.cc
    construct_voronoi_cell.cc
)

target_include_directories(${PROJECT_NAME}-meshing PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(${PROJECT_NAME}-meshing
                        PUBLIC cxx_std_11
)

find_package(Gmsh 4.8.4 REQUIRED MODULE)
find_package(Voro++ REQUIRED MODULE)
find_package(Eigen3 3.3 CONFIG REQUIRED)
find_package(range-v3 0.11.0 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-meshing 
                            PUBLIC external::gmsh
                            PUBLIC external::voro++
                            PUBLIC Eigen3::Eigen
                            PUBLIC range-v3::range-v3)

add_library(${PROJECT_NAME}::meshing ALIAS ${PROJECT_NAME}-meshing)

ae108_install_library(meshing)