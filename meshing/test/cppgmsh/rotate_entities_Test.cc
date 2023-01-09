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
#include "ae108/meshing/cppgmsh/rotate_entities.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct rotate_entities_Test : Test {};

TEST_F(rotate_entities_Test, rotates_rectangle_by_30_degrees) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box = construct_rectangle({0, 0, 0}, {1, 1});
  rotate_entities({box}, {0., 0., 0.}, {0., 0., 1.}, 30. * M_PI / 180.);
  synchronize();

  const auto points = get_points_of(box);
  std::vector<std::array<double, 3>> coords;
  for (const auto &point : points)
    coords.push_back(get_coords_of(point));

  ASSERT_THAT(coords[0], ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  ASSERT_THAT(coords[1][1], DoubleEq(1. / 2.));
  ASSERT_THAT(coords[3][0], DoubleEq(-1. / 2.));
}

TEST_F(rotate_entities_Test, rotates_rectangle_by_pi_sixths) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box = construct_rectangle({0, 0, 0}, {1, 1});
  rotate_entities({box}, {0., 0., 0.}, {0., 0., 1.}, M_PI / 6.);
  synchronize();

  const auto points = get_points_of(box);
  std::vector<std::array<double, 3>> coords;
  for (const auto &point : points)
    coords.push_back(get_coords_of(point));

  ASSERT_THAT(coords[0], ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  ASSERT_THAT(coords[1][1], DoubleEq(1. / 2.));
  ASSERT_THAT(coords[3][0], DoubleEq(-1. / 2.));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108