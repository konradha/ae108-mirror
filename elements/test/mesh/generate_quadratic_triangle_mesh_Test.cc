// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/mesh/generate_quadratic_triangle_mesh.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::IsEmpty;
using testing::SizeIs;
using testing::Test;

namespace ae108 {
namespace elements {
namespace mesh {
namespace {

struct generate_quadratic_triangle_mesh_Test : Test {};

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_mesh_for_zero_granularity) {
  const auto mesh = generate_quadratic_triangle_mesh({{1., 1.}}, {{0, 0}});

  EXPECT_THAT(mesh.connectivity(), SizeIs(0));
  EXPECT_THAT(mesh.number_of_positions(), Eq(0));
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_connectivity_for_granularity_1_1) {
  const auto mesh = generate_quadratic_triangle_mesh({{1., 1.}}, {{1, 1}});

  ASSERT_THAT(mesh.connectivity(), SizeIs(2));

  EXPECT_THAT(mesh.connectivity().at(0), ElementsAre(0, 6, 2, 3, 4, 1));
  EXPECT_THAT(mesh.connectivity().at(1), ElementsAre(6, 8, 2, 7, 5, 4));
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_positions_for_granularity_1_1) {
  const auto mesh = generate_quadratic_triangle_mesh({{1., 1.}}, {{1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(9));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(0.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(.5), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(.5), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(.5), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(1.), DoubleEq(1.)));
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_positions_for_x_size_2) {
  const auto mesh = generate_quadratic_triangle_mesh({{2., 1.}}, {{1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(9));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(0.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(1.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(2.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(2.), DoubleEq(1.)));
  ;
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_positions_for_y_size_2) {
  const auto mesh = generate_quadratic_triangle_mesh({{1., 2.}}, {{1, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(9));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(.5), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(.5), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(.5), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(1.), DoubleEq(2.)));
  ;
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_connectivity_for_granularity_2_1) {
  const auto mesh = generate_quadratic_triangle_mesh({{2., 1.}}, {{2, 1}});

  ASSERT_THAT(mesh.connectivity(), SizeIs(4));
  EXPECT_THAT(mesh.connectivity().at(0), ElementsAre(0, 6, 2, 3, 4, 1));
  EXPECT_THAT(mesh.connectivity().at(1), ElementsAre(6, 8, 2, 7, 5, 4));
  EXPECT_THAT(mesh.connectivity().at(2), ElementsAre(6, 12, 8, 9, 10, 7));
  EXPECT_THAT(mesh.connectivity().at(3), ElementsAre(12, 14, 8, 13, 11, 10));
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_positions_for_granularity_2_1) {
  const auto mesh = generate_quadratic_triangle_mesh({{2., 1.}}, {{2, 1}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(15));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(0.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(.5), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(.5), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(.5), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(1.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(1.5), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.5), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(1.5), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(2.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(2.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(2.), DoubleEq(1.)));
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_connectivity_for_granularity_1_2) {
  const auto mesh = generate_quadratic_triangle_mesh({{1., 2.}}, {{1, 2}});

  ASSERT_THAT(mesh.connectivity(), SizeIs(4));
  EXPECT_THAT(mesh.connectivity().at(0), ElementsAre(0, 10, 2, 5, 6, 1));
  EXPECT_THAT(mesh.connectivity().at(1), ElementsAre(10, 12, 2, 11, 7, 6));
  EXPECT_THAT(mesh.connectivity().at(2), ElementsAre(2, 12, 4, 7, 8, 3));
  EXPECT_THAT(mesh.connectivity().at(3), ElementsAre(12, 14, 4, 13, 9, 8));
}

TEST_F(generate_quadratic_triangle_mesh_Test,
       returns_correct_positions_for_granularity_1_2) {
  const auto mesh = generate_quadratic_triangle_mesh({{1., 2.}}, {{1, 2}});

  ASSERT_THAT(mesh.number_of_positions(), Eq(15));
  EXPECT_THAT(mesh.position_of_vertex(0),
              ElementsAre(DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(1),
              ElementsAre(DoubleEq(0.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(2),
              ElementsAre(DoubleEq(0.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(3),
              ElementsAre(DoubleEq(0.), DoubleEq(1.5)));
  EXPECT_THAT(mesh.position_of_vertex(4),
              ElementsAre(DoubleEq(0.), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(5),
              ElementsAre(DoubleEq(.5), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(6),
              ElementsAre(DoubleEq(.5), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(7),
              ElementsAre(DoubleEq(.5), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(8),
              ElementsAre(DoubleEq(.5), DoubleEq(1.5)));
  EXPECT_THAT(mesh.position_of_vertex(9),
              ElementsAre(DoubleEq(.5), DoubleEq(2.)));
  EXPECT_THAT(mesh.position_of_vertex(10),
              ElementsAre(DoubleEq(1.), DoubleEq(0.)));
  EXPECT_THAT(mesh.position_of_vertex(11),
              ElementsAre(DoubleEq(1.), DoubleEq(.5)));
  EXPECT_THAT(mesh.position_of_vertex(12),
              ElementsAre(DoubleEq(1.), DoubleEq(1.)));
  EXPECT_THAT(mesh.position_of_vertex(13),
              ElementsAre(DoubleEq(1.), DoubleEq(1.5)));
  EXPECT_THAT(mesh.position_of_vertex(14),
              ElementsAre(DoubleEq(1.), DoubleEq(2.)));
}

} // namespace
} // namespace mesh
} // namespace elements
} // namespace ae108