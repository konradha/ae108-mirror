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
#include "ae108/meshing/cppgmsh/get_physical_groups.h"
#include "ae108/meshing/cppgmsh/set_physical_group_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_physical_groups_Test : Test {};

TEST_F(get_physical_groups_Test, gets_all_groups) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle_1 = construct_rectangle({0, 0, 0}, {1, 1});
  const auto rectangle_2 = construct_rectangle({0, 1, 0}, {.5, .5});

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({10., 10., 10.}, {1., 1., 1.});
  const auto box_3 = construct_box({5., 10., 10.}, {1., 1., 1.});
  synchronize();

  const auto group_1 = set_physical_group_of({rectangle_1, rectangle_2});
  const auto group_2 = set_physical_group_of({box_1, box_2});
  const auto group_3 = set_physical_group_of({box_3});

  auto groups = get_physical_groups();
  ASSERT_THAT(groups.size(), Eq(3));
  ASSERT_THAT(groups[0], Eq(group_1));
  ASSERT_THAT(groups[1], Eq(group_2));
  ASSERT_THAT(groups[2], Eq(group_3));
}

TEST_F(get_physical_groups_Test, gets_groups_of_specified_dimension) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle_1 = construct_rectangle({0, 0, 0}, {1, 1});
  const auto rectangle_2 = construct_rectangle({0, 1, 0}, {.5, .5});

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({10., 10., 10.}, {1., 1., 1.});
  const auto box_3 = construct_box({5., 10., 10.}, {1., 1., 1.});
  synchronize();

  const auto group_1 = set_physical_group_of({rectangle_1, rectangle_2});
  const auto group_2 = set_physical_group_of({box_1, box_2});
  const auto group_3 = set_physical_group_of({box_3});

  auto groups_2D = get_physical_groups(2);
  ASSERT_THAT(groups_2D.size(), Eq(1));
  ASSERT_THAT(groups_2D[0], Eq(group_1));

  auto groups_3D = get_physical_groups(3);
  ASSERT_THAT(groups_3D[0], Eq(group_2));
  ASSERT_THAT(groups_3D[1], Eq(group_3));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108