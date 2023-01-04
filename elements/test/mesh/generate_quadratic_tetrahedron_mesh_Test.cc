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

#include "ae108/elements/mesh/generate_quadratic_tetrahedron_mesh.h"
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

struct generate_quadratic_tetrahedron_mesh_Test : Test {};

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       returns_empty_mesh_for_zero_x_granularity) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{1., 1., 1.}}, {{0, 1, 1}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       returns_empty_mesh_for_zero_y_granularity) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{1., 1., 1.}}, {{1, 0, 1}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       returns_empty_mesh_for_zero_z_granularity) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{1., 1., 1.}}, {{1, 1, 0}});

  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
  EXPECT_THAT(mesh.connectivity(), IsEmpty());
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_positions_for_cube) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 2., 2.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(26));
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
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(15),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(16),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(17),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(18),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(19),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(20),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(21),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(22),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(23),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(24),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(25),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(2.)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_connectivity_for_cube) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{1., 1., 1.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 2, 6, 17, 1, 4, 3, 9, 10, 12),
                          ElementsAre(2, 25, 17, 19, 13, 21, 10, 11, 22, 18),
                          ElementsAre(2, 6, 25, 8, 4, 15, 13, 5, 7, 16),
                          ElementsAre(25, 6, 17, 23, 15, 12, 21, 24, 14, 20),
                          ElementsAre(25, 2, 17, 6, 13, 10, 21, 15, 4, 12)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_positions_for_x_scaled_cube) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{4., 2., 2.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(26));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(4.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(4.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(4.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(4.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(4.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(15),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(16),
              ElementsAre(DoubleEq(4.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(17),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(18),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(19),
              ElementsAre(DoubleEq(4.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(20),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(21),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(22),
              ElementsAre(DoubleEq(4.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(23),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(24),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(25),
              ElementsAre(DoubleEq(4.), DoubleEq(2.), DoubleEq(2.)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_positions_for_y_scaled_cube) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 4., 2.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(26));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(0.), DoubleEq(4.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(4.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(2.), DoubleEq(4.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(0.), DoubleEq(4.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(15),
              ElementsAre(DoubleEq(1.), DoubleEq(4.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(16),
              ElementsAre(DoubleEq(2.), DoubleEq(4.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(17),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(18),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(19),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(20),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(21),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(22),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(23),
              ElementsAre(DoubleEq(0.), DoubleEq(4.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(24),
              ElementsAre(DoubleEq(1.), DoubleEq(4.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(25),
              ElementsAre(DoubleEq(2.), DoubleEq(4.), DoubleEq(2.)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_positions_for_z_scaled_cube) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 2., 4.}}, {{1, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(26));
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
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(15),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(16),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(17),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(18),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(19),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(20),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(21),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(22),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(23),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(24),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(25),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(4.)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_connectivity_for_x_scaled_cube) {
  const auto reference =
      generate_quadratic_tetrahedron_mesh({{2., 2., 2.}}, {{1, 1, 1}})
          .connectivity();

  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{4., 2., 2.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(), Eq(reference));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_connectivity_for_y_scaled_cube) {
  const auto reference =
      generate_quadratic_tetrahedron_mesh({{2., 2., 2.}}, {{1, 1, 1}})
          .connectivity();

  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 4., 2.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(), Eq(reference));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_connectivity_for_z_scaled_cube) {
  const auto reference =
      generate_quadratic_tetrahedron_mesh({{2., 2., 2.}}, {{1, 1, 1}})
          .connectivity();

  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 2., 4.}}, {{1, 1, 1}});

  EXPECT_THAT(mesh.connectivity(), Eq(reference));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_positions_for_x_split_cuboid) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{4., 2., 2.}}, {{2, 1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(43));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(3.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(4.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(3.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(4.), DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(3.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(4.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(15),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(16),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(17),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(18),
              ElementsAre(DoubleEq(3.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(19),
              ElementsAre(DoubleEq(4.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(20),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(21),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(22),
              ElementsAre(DoubleEq(4.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(23),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(24),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(25),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(26),
              ElementsAre(DoubleEq(3.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(27),
              ElementsAre(DoubleEq(4.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(28),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(29),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(30),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(31),
              ElementsAre(DoubleEq(3.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(32),
              ElementsAre(DoubleEq(4.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(33),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(34),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(35),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(36),
              ElementsAre(DoubleEq(3.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(37),
              ElementsAre(DoubleEq(4.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(38),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(39),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(40),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(41),
              ElementsAre(DoubleEq(3.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(42),
              ElementsAre(DoubleEq(4.), DoubleEq(2.), DoubleEq(2.)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_positions_for_y_split_cuboid) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 4., 2.}}, {{1, 2, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(43));
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
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(0.), DoubleEq(3.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.), DoubleEq(3.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(2.), DoubleEq(3.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(0.), DoubleEq(4.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(1.), DoubleEq(4.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(2.), DoubleEq(4.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(15),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(16),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(17),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(18),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(19),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(20),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(21),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(22),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(23),
              ElementsAre(DoubleEq(0.), DoubleEq(3.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(24),
              ElementsAre(DoubleEq(2.), DoubleEq(3.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(25),
              ElementsAre(DoubleEq(0.), DoubleEq(4.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(26),
              ElementsAre(DoubleEq(1.), DoubleEq(4.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(27),
              ElementsAre(DoubleEq(2.), DoubleEq(4.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(28),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(29),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(30),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(31),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(32),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(33),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(34),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(35),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(36),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(37),
              ElementsAre(DoubleEq(0.), DoubleEq(3.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(38),
              ElementsAre(DoubleEq(1.), DoubleEq(3.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(39),
              ElementsAre(DoubleEq(2.), DoubleEq(3.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(40),
              ElementsAre(DoubleEq(0.), DoubleEq(4.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(41),
              ElementsAre(DoubleEq(1.), DoubleEq(4.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(42),
              ElementsAre(DoubleEq(2.), DoubleEq(4.), DoubleEq(2.)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_positions_for_z_split_cuboid) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 2., 4.}}, {{1, 1, 2}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(43));
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
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(15),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(16),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(17),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(18),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(19),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(20),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(21),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(22),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(23),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(24),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(25),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(26),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(27),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(28),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(29),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(30),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(31),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(32),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(33),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(3.)));
  EXPECT_THAT(mesh.position_of_vertex(34),
              ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(35),
              ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(36),
              ElementsAre(DoubleEq(2.), DoubleEq(0.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(37),
              ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(38),
              ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(39),
              ElementsAre(DoubleEq(2.), DoubleEq(1.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(40),
              ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(41),
              ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(4.)));
  EXPECT_THAT(mesh.position_of_vertex(42),
              ElementsAre(DoubleEq(2.), DoubleEq(2.), DoubleEq(4.)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_connectivity_for_x_split_cuboid) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{4., 2., 2.}}, {{2, 1, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 2, 10, 28, 1, 6, 5, 15, 16, 20),
                          ElementsAre(2, 40, 28, 30, 21, 34, 16, 17, 35, 29),
                          ElementsAre(2, 10, 40, 12, 6, 24, 21, 7, 11, 25),
                          ElementsAre(40, 10, 28, 38, 24, 20, 34, 39, 23, 33),
                          ElementsAre(40, 2, 28, 10, 21, 16, 34, 24, 6, 20),
                          ElementsAre(2, 4, 14, 32, 3, 9, 8, 18, 19, 22),
                          ElementsAre(14, 40, 32, 42, 26, 36, 22, 27, 41, 37),
                          ElementsAre(2, 40, 14, 12, 21, 26, 8, 7, 25, 13),
                          ElementsAre(2, 32, 40, 30, 18, 36, 21, 17, 31, 35),
                          ElementsAre(32, 40, 14, 2, 36, 26, 22, 18, 21, 8)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_connectivity_for_y_split_cuboid) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 4., 2.}}, {{1, 2, 1}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 2, 6, 28, 1, 4, 3, 15, 16, 18),
                          ElementsAre(2, 36, 28, 30, 19, 32, 16, 17, 33, 29),
                          ElementsAre(2, 6, 36, 8, 4, 21, 19, 5, 7, 22),
                          ElementsAre(36, 6, 28, 34, 21, 18, 32, 35, 20, 31),
                          ElementsAre(36, 2, 28, 6, 19, 16, 32, 21, 4, 18),
                          ElementsAre(6, 8, 14, 36, 7, 11, 10, 21, 22, 24),
                          ElementsAre(14, 40, 36, 42, 26, 38, 24, 27, 41, 39),
                          ElementsAre(6, 40, 14, 12, 23, 26, 10, 9, 25, 13),
                          ElementsAre(6, 36, 40, 34, 21, 38, 23, 20, 35, 37),
                          ElementsAre(36, 40, 14, 6, 38, 26, 24, 21, 23, 10)));
}

TEST_F(generate_quadratic_tetrahedron_mesh_Test,
       has_correct_connectivity_for_z_split_cuboid) {
  const auto mesh =
      generate_quadratic_tetrahedron_mesh({{2., 2., 4.}}, {{1, 1, 2}});

  EXPECT_THAT(mesh.connectivity(),
              ElementsAre(ElementsAre(0, 2, 6, 17, 1, 4, 3, 9, 10, 12),
                          ElementsAre(2, 25, 17, 19, 13, 21, 10, 11, 22, 18),
                          ElementsAre(2, 6, 25, 8, 4, 15, 13, 5, 7, 16),
                          ElementsAre(25, 6, 17, 23, 15, 12, 21, 24, 14, 20),
                          ElementsAre(25, 2, 17, 6, 13, 10, 21, 15, 4, 12),
                          ElementsAre(17, 19, 25, 36, 18, 22, 21, 27, 28, 30),
                          ElementsAre(25, 40, 36, 42, 32, 38, 30, 33, 41, 39),
                          ElementsAre(17, 40, 25, 23, 29, 32, 21, 20, 31, 24),
                          ElementsAre(17, 36, 40, 34, 27, 38, 29, 26, 35, 37),
                          ElementsAre(36, 40, 25, 17, 38, 32, 30, 27, 29, 21)));
}

} // namespace
} // namespace mesh
} // namespace elements
} // namespace ae108