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

#include "ae108/meshing/bounding_box_of.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace meshing {
struct bounding_box_of_Test : Test {};

TEST_F(bounding_box_of_Test, returns_correct_bounding_box_of_single_point) {
  const auto box = bounding_box_of<std::array<double, 2>>({{1, 2}});
  ASSERT_THAT(box.min, ElementsAre(DoubleEq(1), DoubleEq(2)));
  ASSERT_THAT(box.max, ElementsAre(DoubleEq(1), DoubleEq(2)));
};

TEST_F(bounding_box_of_Test, returns_correct_bounding_box_for_multiple_point) {
  const auto box = bounding_box_of<std::array<double, 3>>(
      {{1, 2, 3}, {-1, 2, -3}, {3, 2, 1}});
  ASSERT_THAT(box.min, ElementsAre(DoubleEq(-1), DoubleEq(2), DoubleEq(-3)));
  ASSERT_THAT(box.max, ElementsAre(DoubleEq(3), DoubleEq(2), DoubleEq(3)));
};

} // namespace meshing
} // namespace ae108