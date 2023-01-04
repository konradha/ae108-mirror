// © 2020 ETH Zurich, Mechanics and Materials Lab
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