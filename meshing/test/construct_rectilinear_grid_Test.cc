// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/meshing/construct_rectilinear_grid.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
struct construct_rectilinear_grid_Test : Test {};

TEST_F(construct_rectilinear_grid_Test, returns_correct_1D_grid) {

  std::array<std::vector<double>, 1> instances = {{{1, 2}}};
  const auto [nodes, segments] = construct_rectilinear_grid(instances);

  ASSERT_THAT(nodes.size(), Eq(2));
  ASSERT_THAT(segments.size(), Eq(1));
  ASSERT_THAT(nodes[0], ElementsAre(DoubleEq(1)));
  ASSERT_THAT(nodes[1], ElementsAre(DoubleEq(2)));
  ASSERT_THAT(segments[0], ElementsAre(Eq(0), Eq(1)));
};

TEST_F(construct_rectilinear_grid_Test, returns_correct_2D_grid) {

  std::array<std::vector<double>, 2> instances = {{{1, 2}, {3, 4}}};
  const auto [nodes, segments] = construct_rectilinear_grid(instances);

  ASSERT_THAT(nodes.size(), Eq(4));
  ASSERT_THAT(segments.size(), Eq(4));
  ASSERT_THAT(nodes[0], ElementsAre(DoubleEq(1), DoubleEq(3)));
  ASSERT_THAT(nodes[1], ElementsAre(DoubleEq(1), DoubleEq(4)));
  ASSERT_THAT(nodes[2], ElementsAre(DoubleEq(2), DoubleEq(3)));
  ASSERT_THAT(nodes[3], ElementsAre(DoubleEq(2), DoubleEq(4)));
  ASSERT_THAT(segments[0], ElementsAre(Eq(0), Eq(1)));
  ASSERT_THAT(segments[1], ElementsAre(Eq(2), Eq(3)));
  ASSERT_THAT(segments[2], ElementsAre(Eq(0), Eq(2)));
  ASSERT_THAT(segments[3], ElementsAre(Eq(1), Eq(3)));
};

TEST_F(construct_rectilinear_grid_Test, returns_correct_3D_grid) {

  std::array<std::vector<double>, 3> instances = {{{1, 2}, {3, 4}, {0}}};
  const auto [nodes, segments] = construct_rectilinear_grid(instances);

  ASSERT_THAT(nodes.size(), Eq(4));
  ASSERT_THAT(segments.size(), Eq(4));
  ASSERT_THAT(nodes[0], ElementsAre(DoubleEq(1), DoubleEq(3), DoubleEq(0)));
  ASSERT_THAT(nodes[1], ElementsAre(DoubleEq(1), DoubleEq(4), DoubleEq(0)));
  ASSERT_THAT(nodes[2], ElementsAre(DoubleEq(2), DoubleEq(3), DoubleEq(0)));
  ASSERT_THAT(nodes[3], ElementsAre(DoubleEq(2), DoubleEq(4), DoubleEq(0)));
  ASSERT_THAT(segments[0], ElementsAre(Eq(0), Eq(1)));
  ASSERT_THAT(segments[1], ElementsAre(Eq(2), Eq(3)));
  ASSERT_THAT(segments[2], ElementsAre(Eq(0), Eq(2)));
  ASSERT_THAT(segments[3], ElementsAre(Eq(1), Eq(3)));
};

} // namespace meshing
} // namespace ae108