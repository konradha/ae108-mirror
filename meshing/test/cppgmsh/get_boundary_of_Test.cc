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
#include "ae108/meshing/cppgmsh/get_boundary_of.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_boundary_of_Test : Test {};

TEST_F(get_boundary_of_Test, box_correctly) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box = construct_box({0., 0., 0.}, {1., 1., 1.});

  const auto boundary = get_boundary_of({box});

  ASSERT_THAT(boundary.size(), Eq(6));
  for (std::size_t i = 0; i < boundary.size(); i++) {
    ASSERT_THAT(boundary[i].first, Eq(2));
  }
}

TEST_F(get_boundary_of_Test, rectangle_correctly) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle = construct_rectangle({0., 0., 0.}, {1., 1.});

  const auto boundary = get_boundary_of({rectangle});

  ASSERT_THAT(boundary.size(), Eq(4));
  for (std::size_t i = 0; i < boundary.size(); i++) {
    ASSERT_THAT(boundary[i].first, Eq(1));
  }
}

TEST_F(get_boundary_of_Test, line_correctly) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto point_1 = gmsh::model::occ::addPoint(0., 0., 0.);
  const auto point_2 = gmsh::model::occ::addPoint(1., 1., 1.);

  const auto line = gmsh::model::occ::addLine(point_1, point_2);

  const auto boundary = get_boundary_of({{1, line}});

  ASSERT_THAT(boundary.size(), Eq(2));
  for (std::size_t i = 0; i < boundary.size(); i++) {
    ASSERT_THAT(boundary[i].first, Eq(0));
  }
  ASSERT_THAT(boundary[0].second, Eq(point_1));
  ASSERT_THAT(boundary[1].second, Eq(point_2));
}

TEST_F(get_boundary_of_Test, recursive_rectangle_correctly) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle = construct_rectangle({0., 0., 0.}, {1., 1.});

  const auto boundary = get_boundary_of({rectangle}, true);

  ASSERT_THAT(boundary.size(), Eq(4));
  for (std::size_t i = 0; i < boundary.size(); i++) {
    ASSERT_THAT(boundary[i].first, Eq(0));
  }
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108