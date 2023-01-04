// © 2021 ETH Zurich, Mechanics and Materials Lab
//
// This file is part of ae108.
//
// ae108 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any
// later version.
//
// ae108 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ae108. If not, see <https://www.gnu.org/licenses/>.

#include "ae108/elements/mesh/generate_cuboid_mesh.h"
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

struct generate_cuboid_mesh_Test : Test {};

TEST_F(generate_cuboid_mesh_Test, returns_empty_mesh_for_zero_x_granularity) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 1.}}, {{0, 1, 1}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_cuboid_mesh_Test, returns_empty_mesh_for_zero_y_granularity) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 1.}}, {{1, 0, 1}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_cuboid_mesh_Test, returns_empty_mesh_for_zero_z_granularity) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 1.}}, {{1, 1, 0}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_cuboid_mesh_Test, has_correct_positions_for_cube) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 1.}}, {{1, 1, 1}});

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

TEST_F(generate_cuboid_mesh_Test, has_correct_connectivity_for_cube) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 1.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 3, 2, 4, 5, 7, 6)));
}

TEST_F(generate_cuboid_mesh_Test, has_correct_positions_for_x_scaled_cube) {
  const auto mesh = generate_cuboid_mesh({{2., 1., 1.}}, {{1, 1, 1}});

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

TEST_F(generate_cuboid_mesh_Test, has_correct_positions_for_y_scaled_cube) {
  const auto mesh = generate_cuboid_mesh({{1., 2., 1.}}, {{1, 1, 1}});

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

TEST_F(generate_cuboid_mesh_Test, has_correct_positions_for_z_scaled_cube) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 2.}}, {{1, 1, 1}});

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

TEST_F(generate_cuboid_mesh_Test, has_correct_connectivity_for_x_scaled_cube) {
  const auto mesh = generate_cuboid_mesh({{2., 1., 1.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 3, 2, 4, 5, 7, 6)));
  ;
}

TEST_F(generate_cuboid_mesh_Test, has_correct_connectivity_for_y_scaled_cube) {
  const auto mesh = generate_cuboid_mesh({{1., 2., 1.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 3, 2, 4, 5, 7, 6)));
}

TEST_F(generate_cuboid_mesh_Test, has_correct_connectivity_for_z_scaled_cube) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 2.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 3, 2, 4, 5, 7, 6)));
}

TEST_F(generate_cuboid_mesh_Test, has_correct_positions_for_x_split_cuboid) {
  const auto mesh = generate_cuboid_mesh({{2., 1., 1.}}, {{2, 1, 1}});

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

TEST_F(generate_cuboid_mesh_Test, has_correct_positions_for_y_split_cuboid) {
  const auto mesh = generate_cuboid_mesh({{1., 2., 1.}}, {{1, 2, 1}});

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

TEST_F(generate_cuboid_mesh_Test, has_correct_positions_for_z_split_cuboid) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 2.}}, {{1, 1, 2}});

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

TEST_F(generate_cuboid_mesh_Test, has_correct_connectivity_for_x_split_cuboid) {
  const auto mesh = generate_cuboid_mesh({{2., 1., 1.}}, {{2, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 4, 3, 6, 7, 10, 9),
                          ElementsAre(1, 2, 5, 4, 7, 8, 11, 10)));
}

TEST_F(generate_cuboid_mesh_Test, has_correct_connectivity_for_y_split_cuboid) {
  const auto mesh = generate_cuboid_mesh({{1., 2., 1.}}, {{1, 2, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 3, 2, 6, 7, 9, 8),
                          ElementsAre(2, 3, 5, 4, 8, 9, 11, 10)));
}

TEST_F(generate_cuboid_mesh_Test, has_correct_connectivity_for_z_split_cuboid) {
  const auto mesh = generate_cuboid_mesh({{1., 1., 2.}}, {{1, 1, 2}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 1, 3, 2, 4, 5, 7, 6),
                          ElementsAre(4, 5, 7, 6, 8, 9, 11, 10)));
}

} // namespace
} // namespace mesh
} // namespace elements
} // namespace ae108