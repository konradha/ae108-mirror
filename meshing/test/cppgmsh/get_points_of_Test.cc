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
#include "ae108/meshing/cppgmsh/get_points_of.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct points_of_Test : Test {};

TEST_F(points_of_Test, returns_correct_points_of_point) {
  const auto gmshContext = Context(0, 0);
  const auto tag = gmsh::model::occ::addPoint(0, 0, 0);
  gmsh::model::occ::synchronize();
  const auto points = get_points_of({0, tag});
  ASSERT_THAT(points.size(), Eq(1));
  ASSERT_THAT(points[0], Eq(tag));
}

TEST_F(points_of_Test, returns_correct_points_of_line) {
  const auto gmshContext = Context(0, 0);
  const auto A = gmsh::model::occ::addPoint(0, 0, 0);
  const auto B = gmsh::model::occ::addPoint(1, 0, 0);
  const auto tag = gmsh::model::occ::addLine(A, B);
  gmsh::model::occ::synchronize();
  const auto points = get_points_of({1, tag});
  ASSERT_THAT(points.size(), Eq(2));
  ASSERT_THAT(points[0], Eq(A));
  ASSERT_THAT(points[1], Eq(B));
}

TEST_F(points_of_Test, returns_correct_points_of_rectangle) {
  const auto gmshContext = Context(0, 0);
  const auto tag = gmsh::model::occ::addRectangle(0, 0, 0, 1, 1);
  gmsh::model::occ::synchronize();
  const auto points = get_points_of({2, tag});
  ASSERT_THAT(points.size(), Eq(4));
}

TEST_F(points_of_Test, returns_correct_points_of_box) {
  const auto gmshContext = Context(0, 0);
  const auto tag = gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1);
  gmsh::model::occ::synchronize();
  const auto points = get_points_of({3, tag});
  ASSERT_THAT(points.size(), Eq(8));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108