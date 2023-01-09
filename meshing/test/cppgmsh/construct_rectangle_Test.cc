// © 2022 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include "ae108/meshing/cppgmsh/get_coords_of.h"
#include "ae108/meshing/cppgmsh/get_points_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct construct_rectangle_Test : Test {};

TEST_F(construct_rectangle_Test, constructs_rectangle) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto origin = std::array<double, 3>{0, 0, 0};
  const auto side_lengths = std::array<double, 2>{1, 2};

  const auto rectangle = construct_rectangle(origin, side_lengths);

  ASSERT_THAT(rectangle.first, Eq(2));
  ASSERT_THAT(rectangle.second, Eq(1));

  synchronize();
  const auto points = get_coords_of(get_points_of(rectangle));
  ASSERT_THAT(points[0], ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  ASSERT_THAT(points[1], ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)));
  ASSERT_THAT(points[2], ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(0.)));
  ASSERT_THAT(points[3], ElementsAre(DoubleEq(0.), DoubleEq(2.), DoubleEq(0.)));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108