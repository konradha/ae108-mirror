// © 2021 ETH Zurich, Mechanics and Materials Lab
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
typename Element::Vector create_reference_element_axis() noexcept {
  auto element_axis = typename Element::Vector();
  std::fill(element_axis.begin(), element_axis.end(), 0.);
  element_axis[0] = 1.;
  return element_axis;
}

template <class Element>
typename Element::Vector create_rotated_and_stretched_element_axis() noexcept {
  auto element_axis = typename Element::Vector();
  std::iota(element_axis.begin(), element_axis.end(), 3.);
  return element_axis;
}

template <class ValueType_, std::size_t Dimension_>
inline typename timoshenko::Properties<ValueType_, Dimension_>
create_element_properties() noexcept {
  exit(1);
}

template <>
inline typename timoshenko::Properties<double, 3>
create_element_properties<double, 3>() noexcept {
  return {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
}
template <>
inline typename timoshenko::Properties<double, 2>
create_element_properties<double, 2>() noexcept {
  return {1., 1., 1., 1., 1., 1., 1.};
}

template <class Element_> struct ReferenceConfiguration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(create_reference_element_axis<Element_>(),
                   create_element_properties<typename Element_::value_type,
                                             Element_::dimension()>());
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

template <class Element_> struct RotatedAndStretchedConfiguration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(create_rotated_and_stretched_element_axis<Element_>(),
                   create_element_properties<typename Element_::value_type,
                                             Element_::dimension()>());
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Configurations =
    Types<ReferenceConfiguration<timoshenko::BeamElement<2>>,
          ReferenceConfiguration<timoshenko::BeamElement<3>>,
          RotatedAndStretchedConfiguration<timoshenko::BeamElement<2>>,
          RotatedAndStretchedConfiguration<timoshenko::BeamElement<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElement_Test, Element_Test,
                              Configurations);

struct TimoshenkoBeamElement2D_Test : Test {
  using Element = timoshenko::BeamElement<2>;
  const Element element = ReferenceConfiguration<Element>::create_element();
};

TEST_F(TimoshenkoBeamElement2D_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{1., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(.5));
}

TEST_F(TimoshenkoBeamElement2D_Test,
       computes_energy_with_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{0., 1., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(6. / 13.));
}

TEST_F(TimoshenkoBeamElement2D_Test,
       computes_energy_with_rotational_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{0., 0., 1.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(8. / 13.));
}

struct TimoshenkoBeamElement3D_Test : Test {
  using Element = timoshenko::BeamElement<3>;
  const Element element = ReferenceConfiguration<Element>::create_element();
  const Element::Vector axis = create_reference_element_axis<Element>();
  const Element::value_type L = tensor::as_vector(&axis).norm();
};

TEST_F(TimoshenkoBeamElement3D_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{1., 0., 0., 0., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(.5));
}

TEST_F(TimoshenkoBeamElement3D_Test,
       computes_energy_with_lateral_displacement1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{0., 1., 0., 0., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(6. / 13.));
}

TEST_F(TimoshenkoBeamElement3D_Test,
       computes_energy_with_lateral_displacement2) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{0., 0., 1., 0., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(6. / 13.));
}

TEST_F(TimoshenkoBeamElement3D_Test,
       computes_energy_with_rotational_displacement1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{0., 0., 0., 1., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(1. / 2 / L));
}

TEST_F(TimoshenkoBeamElement3D_Test,
       computes_energy_with_rotational_displacement2) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{0., 0., 0., 0., 1., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq((4. + 12. / L / L) / (1 + 12. / L / L) / L / 2));
}

TEST_F(TimoshenkoBeamElement3D_Test,
       computes_energy_with_rotational_displacement3) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{0., 0., 0., 0., 0., 1.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq((4. + 12. / L / L) / (1 + 12. / L / L) / L / 2));
}

struct TimoshenkoBeamElement2D_rotated_Test : Test {
  using Element = timoshenko::BeamElement<2>;
  const Element element =
      RotatedAndStretchedConfiguration<Element>::create_element();
  const Element::Vector axis =
      create_rotated_and_stretched_element_axis<Element>();
  const Element::value_type L = tensor::as_vector(&axis).norm();
};

TEST_F(TimoshenkoBeamElement2D_rotated_Test,
       computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{3. / L, 4. / L, 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(1. / 2 / L));
}

TEST_F(TimoshenkoBeamElement2D_rotated_Test,
       computes_energy_with_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{4. / L, -3. / L, 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq(6. / (L * L * L + 12 * L)));
}

TEST_F(TimoshenkoBeamElement2D_rotated_Test,
       computes_energy_with_rotational_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{0., 0., 1.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq((4. + 12. / L / L) / (1 + 12. / L / L) / L / 2));
}

struct TimoshenkoBeamElement3D_rotated_Test : Test {
  using Element = timoshenko::BeamElement<3>;
  const Element element =
      RotatedAndStretchedConfiguration<Element>::create_element();
  const Element::Vector axis =
      create_rotated_and_stretched_element_axis<Element>();
  const Element::value_type L = tensor::as_vector(&axis).norm();
};

TEST_F(TimoshenkoBeamElement3D_rotated_Test,
       computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{3. / L, 4. / L, 5. / L, 0., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq(1. / 2. / L));
}

TEST_F(TimoshenkoBeamElement3D_rotated_Test,
       computes_energy_with_lateral_displacement1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{-0.8, 0.6, 0., 0., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq(6. / (L * L * L + 12 * L)));
}

TEST_F(TimoshenkoBeamElement3D_rotated_Test,
       computes_energy_with_rotational_displacement1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{0., 0., 0., 3. / L, 4. / L, 5. / L}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(0.5 / L));
}

TEST_F(TimoshenkoBeamElement3D_rotated_Test,
       computes_energy_with_rotational_displacement2) {
  const auto time = Element::Time{0.};

  double b3 = 1. / sqrt(41. / 16);

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0., 0., 0., 0.}},
      {{0., 0., 0., -0.8, 0.6, 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq((4. + 12. / L / L) / (1 + 12. / L / L) / L / 2));
}

} // namespace
} // namespace elements
} // namespace ae108