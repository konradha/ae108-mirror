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
#include "ae108/meshing/cppgmsh/get_coords_of.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_coords_of_Test : Test {};

TEST_F(get_coords_of_Test, returns_correct_coords_of_point) {
  const auto gmshContext = Context(0, 0);
  const auto tag = gmsh::model::occ::addPoint(1, 2, 3);
  gmsh::model::occ::synchronize();
  const auto coords = get_coords_of(tag);
  ASSERT_THAT(coords[0], Eq(1));
  ASSERT_THAT(coords[1], Eq(2));
  ASSERT_THAT(coords[2], Eq(3));
}

TEST_F(get_coords_of_Test, returns_correct_coords_of_points) {
  const auto gmshContext = Context(0, 0);
  const auto tag_1 = gmsh::model::occ::addPoint(1, 2, 3);
  const auto tag_2 = gmsh::model::occ::addPoint(4, 5, 6);
  gmsh::model::occ::synchronize();
  const auto coords = get_coords_of({tag_1, tag_2});
  ASSERT_THAT(coords[0][0], Eq(1));
  ASSERT_THAT(coords[0][1], Eq(2));
  ASSERT_THAT(coords[0][2], Eq(3));
  ASSERT_THAT(coords[1][0], Eq(4));
  ASSERT_THAT(coords[1][1], Eq(5));
  ASSERT_THAT(coords[1][2], Eq(6));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108