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

#include "ae108/meshing/cppgmsh/Context.h"
#include "ae108/meshing/cppgmsh/construct_box.h"
#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include "ae108/meshing/cppgmsh/set_physical_group_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct set_physical_group_of_Test : Test {};

TEST_F(set_physical_group_of_Test, adds_1D_group) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto line_1 = gmsh::model::occ::addLine(
      gmsh::model::occ::addPoint(0, 0, 0), gmsh::model::occ::addPoint(1, 0, 0));
  const auto line_2 = gmsh::model::occ::addLine(
      gmsh::model::occ::addPoint(1, 0, 0), gmsh::model::occ::addPoint(1, 1, 0));
  synchronize();

  const auto group = set_physical_group_of({{1, line_1}, {1, line_2}});
  ASSERT_THAT(group, Eq(std::pair<int, int>{1, 1}));
}

TEST_F(set_physical_group_of_Test, adds_2D_group) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle_1 = construct_rectangle({0, 0, 0}, {1, 1});
  const auto rectangle_2 = construct_rectangle({0, 1, 0}, {0.5, 0.5});
  synchronize();

  const auto group = set_physical_group_of({rectangle_1, rectangle_2});
  ASSERT_THAT(group, Eq(std::pair<int, int>{2, 1}));
}

TEST_F(set_physical_group_of_Test, adds_3D_group) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({10., 10., 10.}, {1., 1., 1.});
  synchronize();

  const auto group = set_physical_group_of({box_1, box_2});
  ASSERT_THAT(group, Eq(std::pair<int, int>{3, 1}));
}

TEST_F(set_physical_group_of_Test, adds_group_with_one_member) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  synchronize();

  const auto group = set_physical_group_of({box_1});
  ASSERT_THAT(group, Eq(std::pair<int, int>{3, 1}));
}

TEST_F(set_physical_group_of_Test, adds_group_with_three_members) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({10., 10., 10.}, {1., 1., 1.});
  const auto box_3 = construct_box({5., 10., 10.}, {1., 1., 1.});
  synchronize();

  const auto group = set_physical_group_of({box_1, box_2, box_3});
  ASSERT_THAT(group, Eq(std::pair<int, int>{3, 1}));
}
} // namespace cppgmsh
} // namespace meshing
} // namespace ae108