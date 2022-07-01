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
using testing::Eq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace {

template <class Element>
tensor::Tensor<double, Element::dimension()>
create_reference_element_axis() noexcept {
  auto element_axis = tensor::Tensor<double, Element::dimension()>();
  element_axis[0] = 1.;
  return element_axis;
}

template <class Element>
tensor::Tensor<double, Element::dimension()>
create_rotated_and_stretched_element_axis() noexcept {
  auto element_axis = tensor::Tensor<double, Element::dimension()>();
  std::iota(element_axis.begin(), element_axis.end(), 3.);
  return element_axis;
}

template <std::size_t Dimension_>
TimoshenkoBeamProperties<double, Dimension_>
create_element_properties() noexcept;

template <>
TimoshenkoBeamProperties<double, 3> create_element_properties<3>() noexcept {
  return {1., 1., 1., 1., 1., 1., 1., 1.};
}

template <>
TimoshenkoBeamProperties<double, 2> create_element_properties<2>() noexcept {
  return {1., 1., 1., 1., 1.};
}

template <class Element_> struct ReferenceConfiguration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(timoshenko_beam_stiffness_matrix(
        create_reference_element_axis<Element_>(),
        create_element_properties<Element_::dimension()>()));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

template <class Element_> struct RotatedAndStretchedConfiguration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(timoshenko_beam_stiffness_matrix(
        create_rotated_and_stretched_element_axis<Element_>(),
        create_element_properties<Element_::dimension()>()));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Configurations =
    Types<ReferenceConfiguration<TimoshenkoBeamElement<2>>,
          ReferenceConfiguration<TimoshenkoBeamElement<3>>,
          RotatedAndStretchedConfiguration<TimoshenkoBeamElement<2>>,
          RotatedAndStretchedConfiguration<TimoshenkoBeamElement<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElement_Test, Element_Test,
                              Configurations);

struct TimoshenkoBeamElement2D_Test : Test {
  using Element = TimoshenkoBeamElement<2>;
  const Element element = ReferenceConfiguration<Element>::create_element();
};

TEST_F(TimoshenkoBeamElement2D_Test, dimension_is_2) {
  EXPECT_THAT(element.dimension(), Eq(2));
}

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
  using Element = TimoshenkoBeamElement<3>;
  const Element element = ReferenceConfiguration<Element>::create_element();
  const tensor::Tensor<double, 3> axis =
      create_reference_element_axis<Element>();
  const double L = tensor::as_vector(&axis).norm();
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

TEST_F(TimoshenkoBeamElement3D_Test, dimension_is_3) {
  EXPECT_THAT(element.dimension(), Eq(3));
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
  using Element = TimoshenkoBeamElement<2>;
  const Element element =
      RotatedAndStretchedConfiguration<Element>::create_element();
  const tensor::Tensor<double, 2> axis =
      create_rotated_and_stretched_element_axis<Element>();
  const double L = tensor::as_vector(&axis).norm();
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
  using Element = TimoshenkoBeamElement<3>;
  const Element element =
      RotatedAndStretchedConfiguration<Element>::create_element();
  const tensor::Tensor<double, 3> axis =
      create_rotated_and_stretched_element_axis<Element>();
  const double L = tensor::as_vector(&axis).norm();
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