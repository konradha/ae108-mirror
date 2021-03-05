// © 2021 ETH Zurich, Mechanics and Materials Lab
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "ae108/elements/mesh/generate_tetrahedron_mesh.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::IsEmpty;
using testing::Test;

namespace ae108 {
namespace elements {
namespace mesh {
namespace {

struct generate_tetrahedron_mesh_Test : Test {};

TEST_F(generate_tetrahedron_mesh_Test,
       returns_empty_mesh_for_zero_x_granularity) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 1.}}, {{0, 1, 1}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_tetrahedron_mesh_Test,
       returns_empty_mesh_for_zero_y_granularity) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 1.}}, {{1, 0, 1}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_tetrahedron_mesh_Test,
       returns_empty_mesh_for_zero_z_granularity) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 1.}}, {{1, 1, 0}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_tetrahedron_mesh_Test, has_correct_positions_for_cube) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 1.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(8));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(1.)));
}

TEST_F(generate_tetrahedron_mesh_Test, has_correct_connectivity_for_cube) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 1.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 2, 4), ElementsAre(1, 7, 4, 5),
                          ElementsAre(1, 2, 7, 3), ElementsAre(7, 2, 4, 6),
                          ElementsAre(7, 1, 4, 2)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_positions_for_x_scaled_cube) {
  const auto mesh = generate_tetrahedron_mesh({{2., 1., 1.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(8));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(1.)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_positions_for_y_scaled_cube) {
  const auto mesh = generate_tetrahedron_mesh({{1., 2., 1.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(8));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(1.)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_positions_for_z_scaled_cube) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 2.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(8));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(2.)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_connectivity_for_x_scaled_cube) {
  const auto mesh = generate_tetrahedron_mesh({{2., 1., 1.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 2, 4), ElementsAre(1, 7, 4, 5),
                          ElementsAre(1, 2, 7, 3), ElementsAre(7, 2, 4, 6),
                          ElementsAre(7, 1, 4, 2)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_connectivity_for_y_scaled_cube) {
  const auto mesh = generate_tetrahedron_mesh({{1., 2., 1.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 2, 4), ElementsAre(1, 7, 4, 5),
                          ElementsAre(1, 2, 7, 3), ElementsAre(7, 2, 4, 6),
                          ElementsAre(7, 1, 4, 2)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_connectivity_for_z_scaled_cube) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 2.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 2, 4), ElementsAre(1, 7, 4, 5),
                          ElementsAre(1, 2, 7, 3), ElementsAre(7, 2, 4, 6),
                          ElementsAre(7, 1, 4, 2)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_positions_for_x_split_cuboid) {
  const auto mesh = generate_tetrahedron_mesh({{2., 1., 1.}}, {{2, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(12));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(1.)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_positions_for_y_split_cuboid) {
  const auto mesh = generate_tetrahedron_mesh({{1., 2., 1.}}, {{1, 2, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(12));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(1.)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_positions_for_z_split_cuboid) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 2.}}, {{1, 1, 2}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(12));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(2.)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_connectivity_for_x_split_cuboid) {
  const auto mesh = generate_tetrahedron_mesh({{2., 1., 1.}}, {{2, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 3, 6), ElementsAre(1, 10, 6, 7),
                          ElementsAre(1, 3, 10, 4), ElementsAre(10, 3, 6, 9),
                          ElementsAre(10, 1, 6, 3), ElementsAre(1, 2, 5, 8),
                          ElementsAre(5, 10, 8, 11), ElementsAre(1, 10, 5, 4),
                          ElementsAre(1, 8, 10, 7), ElementsAre(8, 10, 5, 1)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_connectivity_for_y_split_cuboid) {
  const auto mesh = generate_tetrahedron_mesh({{1., 2., 1.}}, {{1, 2, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 2, 6), ElementsAre(1, 9, 6, 7),
                          ElementsAre(1, 2, 9, 3), ElementsAre(9, 2, 6, 8),
                          ElementsAre(9, 1, 6, 2), ElementsAre(2, 3, 5, 9),
                          ElementsAre(5, 10, 9, 11), ElementsAre(2, 10, 5, 4),
                          ElementsAre(2, 9, 10, 8), ElementsAre(9, 10, 5, 2)));
}

TEST_F(generate_tetrahedron_mesh_Test,
       has_correct_connectivity_for_z_split_cuboid) {
  const auto mesh = generate_tetrahedron_mesh({{1., 1., 2.}}, {{1, 1, 2}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 2, 4), ElementsAre(1, 7, 4, 5),
                          ElementsAre(1, 2, 7, 3), ElementsAre(7, 2, 4, 6),
                          ElementsAre(7, 1, 4, 2), ElementsAre(4, 5, 7, 9),
                          ElementsAre(7, 10, 9, 11), ElementsAre(4, 10, 7, 6),
                          ElementsAre(4, 9, 10, 8), ElementsAre(9, 10, 7, 4)));
}

} // namespace
} // namespace mesh
} // namespace elements
} // namespace ae108