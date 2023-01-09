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
#include "ae108/meshing/cppgmsh/get_entities_of.h"
#include "ae108/meshing/cppgmsh/get_points_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include "ae108/meshing/cppgmsh/translate_entities.h"
#include <gmock/gmock.h>
#include <math.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct translate_entities_Test : Test {};

TEST_F(translate_entities_Test, translates_rectangle) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle = construct_rectangle({0, 0, 0}, {1, 1});

  translate_entities({rectangle}, {0.1, 0.2, 0.3});

  synchronize();
  const auto points = get_points_of(rectangle);
  std::vector<std::array<double, 3>> coords;
  for (const auto &point : points)
    coords.push_back(get_coords_of(point));

  ASSERT_THAT(coords[0],
              ElementsAre(DoubleEq(0.1), DoubleEq(0.2), DoubleEq(0.3)));
  ASSERT_THAT(coords[1],
              ElementsAre(DoubleEq(1.1), DoubleEq(0.2), DoubleEq(0.3)));
  ASSERT_THAT(coords[2],
              ElementsAre(DoubleEq(1.1), DoubleEq(1.2), DoubleEq(0.3)));
  ASSERT_THAT(coords[3],
              ElementsAre(DoubleEq(0.1), DoubleEq(1.2), DoubleEq(0.3)));
}

TEST_F(translate_entities_Test, makes_translated_copy_of_rectangle) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle = construct_rectangle({0, 0, 0}, {1, 1});

  const auto translated_rectangle =
      translate_entities({rectangle}, {0.1, 0.2, 0.3}, true)[0];

  synchronize();

  const auto rectangles = get_entities_of(2);
  ASSERT_THAT(rectangles.size(), Eq(2));

  const auto points = get_points_of(translated_rectangle);
  std::vector<std::array<double, 3>> coords;
  for (const auto &point : points)
    coords.push_back(get_coords_of(point));

  ASSERT_THAT(coords[0],
              ElementsAre(DoubleEq(0.1), DoubleEq(0.2), DoubleEq(0.3)));
  ASSERT_THAT(coords[1],
              ElementsAre(DoubleEq(1.1), DoubleEq(0.2), DoubleEq(0.3)));
  ASSERT_THAT(coords[2],
              ElementsAre(DoubleEq(1.1), DoubleEq(1.2), DoubleEq(0.3)));
  ASSERT_THAT(coords[3],
              ElementsAre(DoubleEq(0.1), DoubleEq(1.2), DoubleEq(0.3)));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108