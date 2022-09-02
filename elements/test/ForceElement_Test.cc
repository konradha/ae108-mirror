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
#include "ae108/elements/ForceElement.h"
#include <gmock/gmock.h>
#include <numeric>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace {

/**
 * @brief Returns a force of [1., 2., 3., ...] (depending on the number of
 * degrees of freedom).
 */
template <class Element> typename Element::Force create_test_force() noexcept {
  auto force = typename Element::Force();
  std::iota(force.begin(), force.end(), 1.);
  return force;
}

template <class Element_> struct Configuration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(create_test_force<Element_>());
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Configurations =
    Types<Configuration<ForceElement<1>>, Configuration<ForceElement<2>>>;
INSTANTIATE_TYPED_TEST_CASE_P(ForceElement_Test, Element_Test, Configurations);

struct ForceElement_Test : Test {
  static constexpr auto dimension = 2;

  using Element = ForceElement<dimension>;
  const Element element = Element(create_test_force<Element>());
};

TEST_F(ForceElement_Test, dimension_is_2) {
  EXPECT_THAT(element.dimension(), Eq(2));
}

TEST_F(ForceElement_Test, number_of_degrees_of_freedom_is_2) {
  EXPECT_THAT(element.degrees_of_freedom(), Eq(2));
}

TEST_F(ForceElement_Test, computes_correct_energy) {
  const auto time = Element::Time{0.};
  const auto displacements = Element::NodalDisplacements({{{{7., 77.}}}});

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq(7. * 1. + 77. * 2.));
}

TEST_F(ForceElement_Test, returns_provided_force) {
  const auto time = Element::Time{0.};
  const auto displacements = Element::NodalDisplacements();

  EXPECT_THAT(element.computeForces(displacements, time),
              ElementsAre(ElementsAre(1., 2.)));
}

} // namespace
} // namespace elements
} // namespace ae108