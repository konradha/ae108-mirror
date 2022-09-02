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

#include "Element_Test.h"
#include "ae108/elements/Minimal.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Eq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace {

template <class Element_> struct Configuration {
  using Element = Element_;
  static Element create_element() noexcept { return Element(); }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Configurations =
    Types<Configuration<Minimal<1, 1>>, Configuration<Minimal<2, 1>>,
          Configuration<Minimal<2, 2>>>;
INSTANTIATE_TYPED_TEST_CASE_P(Minimal_Test, Element_Test, Configurations);

struct Minimal_Test : Test {
  static constexpr auto size = 2;
  static constexpr auto dimension = 1;

  using Element = Minimal<size, dimension>;
  const Element element = Element();
};

TEST_F(Minimal_Test, dimension_is_1) {
  EXPECT_THAT(element.dimension(), Eq(1));
}

TEST_F(Minimal_Test, number_of_degrees_of_freedom_is_1) {
  EXPECT_THAT(element.dimension(), Eq(1));
}

TEST_F(Minimal_Test, computes_zero_energy) {
  const auto time = Element::Time{0.};
  const auto displacements = Element::NodalDisplacements();

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(0.));
}

} // namespace
} // namespace elements
} // namespace ae108