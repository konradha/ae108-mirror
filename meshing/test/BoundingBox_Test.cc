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

#include "ae108/meshing/BoundingBox.h"
#include <array>
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace meshing {
struct BoundingBox_Test : Test {};

TEST_F(BoundingBox_Test, returns_correct_bbox) {
  const auto bbox = BoundingBox<std::array<double, 3>>({{0, 0, 0}, {1, 1, 1}});
  EXPECT_THAT(bbox.min, ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)));
  EXPECT_THAT(bbox.max, ElementsAre(DoubleEq(1.), DoubleEq(1.), DoubleEq(1.)));
};

} // namespace meshing
} // namespace ae108