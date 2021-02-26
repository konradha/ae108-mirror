// © 2020 ETH Zurich, Mechanics and Materials Lab
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "Shape_Test.h"
#include "ae108/elements/shape/Tri3.h"
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

using Configurations = Types<Tri3>;

INSTANTIATE_TYPED_TEST_CASE_P(Tri3_Test, Shape_Test, Configurations);

struct Tri3_Test : Test {
  using Shape = Tri3;
};

TEST_F(Tri3_Test, points_are_correct) {
  const auto points = get_points<Shape>();

  EXPECT_THAT(points, ElementsAre(ElementsAre(DoubleEq(0.), DoubleEq(0.)),
                                  ElementsAre(DoubleEq(1.), DoubleEq(0.)),
                                  ElementsAre(DoubleEq(0.), DoubleEq(1.))));
}

} // namespace
} // namespace shape
} // namespace elements
} // namespace ae108