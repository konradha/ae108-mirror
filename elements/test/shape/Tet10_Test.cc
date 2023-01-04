// © 2021 ETH Zurich, Mechanics and Materials Lab
//
// This file is part of ae108.
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

#include "Shape_Test.h"
#include "ae108/elements/shape/Tet10.h"
#include "ae108/elements/shape/get_points.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace shape {
namespace {

using Configurations = Types<Tet10>;

INSTANTIATE_TYPED_TEST_CASE_P(Tet10_Test, Shape_Test, Configurations);

struct Tet10_Test : Test {
  using Shape = Tet10;
};

TEST_F(Tet10_Test, points_are_correct) {
  const auto points = get_points<Shape>();

  EXPECT_THAT(
      points,
      ElementsAre(ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(0.)),
                  ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.)),
                  ElementsAre(DoubleEq(0.), DoubleEq(1.), DoubleEq(0.)),
                  ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(1.)),
                  ElementsAre(DoubleEq(.5), DoubleEq(0.), DoubleEq(0.)),
                  ElementsAre(DoubleEq(.5), DoubleEq(.5), DoubleEq(0.)),
                  ElementsAre(DoubleEq(0.), DoubleEq(.5), DoubleEq(0.)),
                  ElementsAre(DoubleEq(0.), DoubleEq(0.), DoubleEq(.5)),
                  ElementsAre(DoubleEq(.5), DoubleEq(0.), DoubleEq(.5)),
                  ElementsAre(DoubleEq(0.), DoubleEq(.5), DoubleEq(.5))));
}

} // namespace
} // namespace shape
} // namespace elements
} // namespace ae108