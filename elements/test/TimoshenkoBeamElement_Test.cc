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
#include "ae108/elements/TimoshenkoBeamElement.h"
#include <gmock/gmock.h>
#include <numeric>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace {

template <class Element>
typename Element::Vector create_node_positions() noexcept {
  auto node_positions = typename Element::Vector();
  std::iota(node_positions.begin(), node_positions.end(), 1.);
  return node_positions;
}

template <class ValueType_, std::size_t Dimension_>
inline typename timoshenko::Properties<ValueType_, Dimension_>
create_element_properties() noexcept {
  exit(1);
}

template <>
inline typename timoshenko::Properties<double, 3>
create_element_properties<double, 3>() noexcept {
  return typename timoshenko::Properties<double, 3>(1., 1., 1., 1., 1., 1., 1.,
                                                    1., 1.);
}
template <>
inline typename timoshenko::Properties<double, 2>
create_element_properties<double, 2>() noexcept {
  return typename timoshenko::Properties<double, 2>(1., 1., 1., 1., 1., 1.);
}

template <class Element_> struct Configuration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(create_node_positions<Element_>(),
                   create_element_properties<typename Element_::value_type,
                                             Element_::dimension()>());
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Configurations = Types<Configuration<timoshenko::BeamElement<2>>,
                             Configuration<timoshenko::BeamElement<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElement_Test, Element_Test,
                              Configurations);

} // namespace
} // namespace elements
} // namespace ae108