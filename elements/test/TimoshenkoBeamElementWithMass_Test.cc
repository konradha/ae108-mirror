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

#include "ae108/elements/TimoshenkoBeamElementWithMass.h"
#include <Eigen/Dense>
#include <gmock/gmock.h>
#include <numeric>

using testing::DoubleEq;
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
  return {1.42, 1.64, 1.46, 1.22, 1.67, 1.98, 1.45, 1.12, 1.97};
}

template <>
TimoshenkoBeamProperties<double, 2> create_element_properties<2>() noexcept {
  return {1.12, 1.75, 1.34, 1.78, 1.06, 1.23};
}

template <class Element_> struct ReferenceConfigurationLumped {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(Element::StiffnessMatrix::Zero(),
                   timoshenko_beam_lumped_mass_matrix(
                       create_reference_element_axis<Element_>(),
                       create_element_properties<Element_::dimension()>()));
  }

  static typename Element::Time calculate_mass() noexcept {
    const auto element_axis = create_reference_element_axis<Element_>();
    const auto properties = create_element_properties<Element_::dimension()>();
    return tensor::as_vector(&element_axis).norm() * properties.area *
           properties.density;
  }
};

template <class Element_> struct ReferenceConfigurationConsistent {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(Element::StiffnessMatrix::Zero(),
                   timoshenko_beam_consistent_mass_matrix(
                       create_reference_element_axis<Element_>(),
                       create_element_properties<Element_::dimension()>()));
  }

  static typename Element::Time calculate_mass() noexcept {
    const auto element_axis = create_reference_element_axis<Element_>();
    const auto properties = create_element_properties<Element_::dimension()>();
    return tensor::as_vector(&element_axis).norm() * properties.area *
           properties.density;
  }
};

template <class Element_> struct RotatedAndStretchedConfigurationConsistent {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(Element::StiffnessMatrix::Zero(),
                   timoshenko_beam_consistent_mass_matrix(
                       create_rotated_and_stretched_element_axis<Element_>(),
                       create_element_properties<Element_::dimension()>()));
  }
  static typename Element::Time calculate_mass() noexcept {
    const auto element_axis =
        create_rotated_and_stretched_element_axis<Element_>();
    const auto properties = create_element_properties<Element_::dimension()>();
    return tensor::as_vector(&element_axis).norm() * properties.area *
           properties.density;
  }
};

template <class Element_> struct RotatedAndStretchedConfigurationLumped {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(Element::StiffnessMatrix::Zero(),
                   timoshenko_beam_lumped_mass_matrix(
                       create_rotated_and_stretched_element_axis<Element_>(),
                       create_element_properties<Element_::dimension()>()));
  }
  static typename Element::Time calculate_mass() noexcept {
    const auto element_axis =
        create_rotated_and_stretched_element_axis<Element_>();
    const auto properties = create_element_properties<Element_::dimension()>();
    return tensor::as_vector(&element_axis).norm() * properties.area *
           properties.density;
  }
};

template <typename TestConfiguration>
struct ElementWithMass_Test : ::testing::Test {
  using Element = typename TestConfiguration::Element;

  Element element = TestConfiguration::create_element();
  const double mass = TestConfiguration::calculate_mass();

  /**
   * Checks that the mass matrix is positive (semi-)definite
   */
  void check_semi_positive_definite() const noexcept {
    const auto mass_matrix = element.computeMassMatrix();
    Eigen::LLT<typename Element::MassMatrix> llt_of_mass_matrix(mass_matrix);

    ASSERT_THAT((mass_matrix.isApprox(mass_matrix.transpose()) &&
                 llt_of_mass_matrix.info() != Eigen::NumericalIssue),
                true);
  }

  /**
   * Checks that the mass matrix complies with Newton's second law
   */
  void check_second_law_of_motion() const noexcept {

    const Eigen::Matrix<double, Element::dimension(), 1>
        rigid_body_translational_acceleration =
            Eigen::Matrix<double, Element::dimension(), 1>::Random();

    Eigen::Matrix<double, Element::size() * Element::degrees_of_freedom(), 1>
        acceleration = Eigen::Matrix<
            double, Element::size() * Element::degrees_of_freedom(), 1>::Zero();
    for (std::size_t node = 0; node < Element::size(); node++)
      for (std::size_t dim = 0; dim < Element::dimension(); dim++)
        acceleration(node * Element::degrees_of_freedom() + dim) +=
            rigid_body_translational_acceleration(dim);

    const auto mass_matrix = element.computeMassMatrix();
    const auto force = mass_matrix * acceleration;

    Eigen::Matrix<double, Element::dimension(), 1> total_force =
        Eigen::Matrix<double, Element::dimension(), 1>::Zero();
    for (std::size_t node = 0; node < Element::size(); node++)
      for (std::size_t dim = 0; dim < Element::dimension(); dim++)
        total_force(dim) += force(node * Element::degrees_of_freedom() + dim);

    EXPECT_THAT(
        total_force.isApprox(mass * rigid_body_translational_acceleration),
        true);
  }
};

TYPED_TEST_CASE_P(ElementWithMass_Test);
TYPED_TEST_P(ElementWithMass_Test, mass_matrix_is_semi_positive_definite) {
  this->check_semi_positive_definite();
}
TYPED_TEST_P(ElementWithMass_Test,
             mass_matrix_complies_with_newtons_second_law) {
  this->check_second_law_of_motion();
}
REGISTER_TYPED_TEST_CASE_P(ElementWithMass_Test,
                           mass_matrix_is_semi_positive_definite,
                           mass_matrix_complies_with_newtons_second_law);

using MassConfigurations = Types<
    ReferenceConfigurationLumped<TimoshenkoBeamElementWithMass<2>>,
    ReferenceConfigurationLumped<TimoshenkoBeamElementWithMass<3>>,
    RotatedAndStretchedConfigurationLumped<TimoshenkoBeamElementWithMass<2>>,
    RotatedAndStretchedConfigurationLumped<TimoshenkoBeamElementWithMass<3>>,
    ReferenceConfigurationConsistent<TimoshenkoBeamElementWithMass<2>>,
    ReferenceConfigurationConsistent<TimoshenkoBeamElementWithMass<3>>,
    RotatedAndStretchedConfigurationConsistent<
        TimoshenkoBeamElementWithMass<2>>,
    RotatedAndStretchedConfigurationConsistent<
        TimoshenkoBeamElementWithMass<3>>>;

INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElementWithMass_Test,
                              ElementWithMass_Test, MassConfigurations);

} // namespace
} // namespace elements
} // namespace ae108