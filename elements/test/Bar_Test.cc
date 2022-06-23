// © 2021, 2022 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/elements/Bar.h"
#include <gmock/gmock.h>
#include <numeric>

using testing::DoubleEq;
using testing::DoubleNear;
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

constexpr static BarProperties<double> properties{.3, .4};

template <class Element_> struct ReferenceConfiguration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(bar_stiffness_matrix<Element::dimension()>(
        create_reference_element_axis<Element_>(), properties));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

template <class Element_> struct RotatedAndStretchedConfiguration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(bar_stiffness_matrix<Element::dimension()>(
        create_rotated_and_stretched_element_axis<Element_>(), properties));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Configurations =
    Types<ReferenceConfiguration<Bar<1>>, ReferenceConfiguration<Bar<2>>,
          ReferenceConfiguration<Bar<3>>,
          RotatedAndStretchedConfiguration<Bar<1>>,
          RotatedAndStretchedConfiguration<Bar<2>>,
          RotatedAndStretchedConfiguration<Bar<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(Bar_Test, Element_Test, Configurations);

struct Bar1D_Test : Test {
  using Element = Bar<1>;
  const Element element = ReferenceConfiguration<Element>::create_element();
};

TEST_F(Bar1D_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0.}},
      {{1.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus));
}

TEST_F(Bar1D_Test, computes_energy_with_double_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0.}},
      {{2.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * 4. * properties.cross_section * properties.young_modulus));
}

struct Bar1D_rotated_Test : Test {
  using Element = Bar<1>;
  const Element element =
      RotatedAndStretchedConfiguration<Element>::create_element();
  const tensor::Tensor<double, 1> axis =
      create_rotated_and_stretched_element_axis<Element>();
  const double L = tensor::as_vector(&axis).norm();
};

TEST_F(Bar1D_rotated_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0.}},
      {{1.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus / L));
}

TEST_F(Bar1D_rotated_Test, computes_energy_with_double_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0.}},
      {{2.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq(.5 * 4. * properties.cross_section *
                       properties.young_modulus / L));
}

struct Bar2D_Test : Test {
  using Element = Bar<2>;
  const Element element = ReferenceConfiguration<Element>::create_element();
};

TEST_F(Bar2D_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{1., 0.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus));
}

TEST_F(Bar2D_Test, computes_energy_with_double_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{2., 0.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * 4. * properties.cross_section * properties.young_modulus));
}

TEST_F(Bar2D_Test, computes_energy_with_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{0., 1.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(0.));
}

TEST_F(Bar2D_Test,
       computes_energy_with_combined_axial_and_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{1., 1.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus));
}

struct Bar3D_Test : Test {
  using Element = Bar<3>;
  const Element element = ReferenceConfiguration<Element>::create_element();
  const tensor::Tensor<double, 3> axis =
      create_reference_element_axis<Element>();
  const double L = tensor::as_vector(&axis).norm();
};

TEST_F(Bar3D_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{axis[0] / L, axis[1] / L, axis[2] / L}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus));
}

TEST_F(Bar3D_Test, computes_energy_with_double_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{2., 0., 0.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * 4. * properties.cross_section * properties.young_modulus));
}

TEST_F(Bar3D_Test, computes_energy_with_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{0., 1., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(0.));
}

TEST_F(Bar3D_Test,
       computes_energy_with_combined_axial_and_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{1., 1., 0.}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus));
}

struct Bar2D_rotated_Test : Test {
  using Element = Bar<2>;
  const Element element =
      RotatedAndStretchedConfiguration<Element>::create_element();
  const tensor::Tensor<double, 2> axis =
      create_rotated_and_stretched_element_axis<Element>();
  const double L = tensor::as_vector(&axis).norm();
};

TEST_F(Bar2D_rotated_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{axis[0] / L, axis[1] / L}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus / L));
}

TEST_F(Bar2D_rotated_Test, computes_energy_with_double_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{2. * axis[0] / L, 2. * axis[1] / L}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq(.5 * 4. * properties.cross_section *
                       properties.young_modulus / L));
}

TEST_F(Bar2D_rotated_Test, computes_energy_with_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{axis[1] / L, -axis[0] / L}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleNear(0., 1e-18));
}

TEST_F(Bar2D_rotated_Test,
       computes_energy_with_combined_axial_and_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{(axis[0] + axis[1]) / L, (axis[1] - axis[0]) / L}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus / L));
}

struct Bar3D_rotated_Test : Test {
  using Element = Bar<3>;
  const Element element =
      RotatedAndStretchedConfiguration<Element>::create_element();
  const tensor::Tensor<double, 3> axis =
      create_rotated_and_stretched_element_axis<Element>();
  const double L = tensor::as_vector(&axis).norm();
};

TEST_F(Bar3D_rotated_Test, computes_energy_with_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{axis[0] / L, axis[1] / L, axis[2] / L}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus / L));
}

TEST_F(Bar3D_rotated_Test, computes_energy_with_double_axial_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{2. * axis[0] / L, 2. * axis[1] / L, 2. * axis[2] / L}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleEq(.5 * 4. * properties.cross_section *
                       properties.young_modulus / L));
}

TEST_F(Bar3D_rotated_Test, computes_energy_with_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{-axis[1] / L, axis[0] / L, 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time),
              DoubleNear(0., 1e-18));
}

TEST_F(Bar3D_rotated_Test,
       computes_energy_with_combined_axial_and_lateral_displacement) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{(axis[0] - axis[1]) / L, (axis[1] + axis[0]) / L, axis[2] / L}},
  }};

  EXPECT_THAT(
      element.computeEnergy(displacements, time),
      DoubleEq(.5 * properties.cross_section * properties.young_modulus / L));
}

} // namespace
} // namespace elements
} // namespace ae108