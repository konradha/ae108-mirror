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
#include "ae108/elements/ElementWithMass.h"
#include "ae108/elements/TimoshenkoBeamElement.h"
#include "ae108/elements/compute_mass_matrix.h"
#include <Eigen/Dense>
#include <gmock/gmock.h>
#include <numeric>

using testing::DoubleEq;
using testing::Eq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace {

template <std::size_t Dimension>
tensor::Tensor<double, Dimension> create_reference_element_axis() noexcept {
  auto element_axis = tensor::Tensor<double, Dimension>();
  element_axis[0] = 1.;
  return element_axis;
}

template <std::size_t Dimension>
tensor::Tensor<double, Dimension>
create_rotated_and_stretched_element_axis() noexcept {
  auto element_axis = tensor::Tensor<double, Dimension>();
  std::iota(element_axis.begin(), element_axis.end(), 3.);
  return element_axis;
}

template <std::size_t Dimension_>
TimoshenkoBeamProperties<double, Dimension_>
create_element_properties() noexcept {
  static_assert(Dimension_ == 2 || Dimension_ == 3,
                "Only dimensions 2 and 3 are supported.");
  if constexpr (Dimension_ == 2) {
    return {1., 1., 1., 1., 1.};
  } else {
    return {1., 1., 1., 1., 1., 1., 1., 1.};
  }
}

template <std::size_t Dimension>
using create_axis_t = decltype(&create_reference_element_axis<Dimension>);

template <class Element_, create_axis_t<Element_::dimension()> CreateAxis_>
struct Configuration {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(timoshenko_beam_stiffness_matrix(
        CreateAxis_(), create_element_properties<Element_::dimension()>()));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Configurations = Types<
    Configuration<TimoshenkoBeamElement<2>, &create_reference_element_axis<2>>,
    Configuration<TimoshenkoBeamElement<3>, &create_reference_element_axis<3>>,
    Configuration<TimoshenkoBeamElement<2>,
                  &create_rotated_and_stretched_element_axis<2>>,
    Configuration<TimoshenkoBeamElement<3>,
                  &create_rotated_and_stretched_element_axis<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElement_Test, Element_Test,
                              Configurations);

struct TimoshenkoBeamElement2D_Test : Test {
  using Element = TimoshenkoBeamElement<2>;
  const Element element =
      Configuration<Element,
                    &create_reference_element_axis<2>>::create_element();
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
  const Element element =
      Configuration<Element,
                    &create_reference_element_axis<3>>::create_element();
  const tensor::Tensor<double, 3> axis =
      create_reference_element_axis<Element::dimension()>();
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
  const Element element = Configuration<
      Element, &create_rotated_and_stretched_element_axis<2>>::create_element();
  const tensor::Tensor<double, 2> axis =
      create_rotated_and_stretched_element_axis<Element::dimension()>();
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
  const Element element = Configuration<
      Element, &create_rotated_and_stretched_element_axis<3>>::create_element();
  const tensor::Tensor<double, 3> axis =
      create_rotated_and_stretched_element_axis<Element::dimension()>();
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

template <std::size_t Dimension_> double create_element_density() noexcept {
  return .123;
}

template <std::size_t Dimension>
using compute_mass_t = decltype(&timoshenko_beam_lumped_mass_matrix<Dimension>);

template <class Element_, create_axis_t<Element_::dimension()> CreateAxis_,
          compute_mass_t<Element_::dimension()> ComputeMass_>
struct ConfigurationWithMass {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(
        typename Element::Element(Element::StiffnessMatrix::Zero()),
        ComputeMass_(CreateAxis_(),
                     create_element_properties<Element_::dimension()>(),
                     create_element_density<Element_::dimension()>()));
  }
  static typename Element::Time calculate_mass() noexcept {
    const auto element_axis = CreateAxis_();
    const auto properties = create_element_properties<Element_::dimension()>();
    const auto density = create_element_density<Element_::dimension()>();
    return tensor::as_vector(&element_axis).norm() * properties.area * density;
  }
};

template <typename TestConfiguration>
struct TimoshenkoBeamElementWithMass_Test : ::testing::Test {
  using Element = typename TestConfiguration::Element;

  const Element element = TestConfiguration::create_element();
  const double mass = TestConfiguration::calculate_mass();

  /**
   * Checks that the mass matrix is positive (semi-)definite.
   */
  void check_semi_positive_definite() const noexcept {
    const auto mass_matrix = compute_mass_matrix(element);
    ASSERT_THAT(mass_matrix.isApprox(mass_matrix.transpose()), Eq(true));

    Eigen::LDLT<typename Element::MassMatrix> ldlt_of_mass_matrix(mass_matrix);
    EXPECT_THAT(ldlt_of_mass_matrix.info() != Eigen::NumericalIssue, Eq(true));
  }

  /**
   * Checks that the mass matrix complies with Newton's second law.
   */
  void check_second_law_of_motion() const noexcept {
    const auto rigid_body_translational_acceleration =
        Eigen::Matrix<double, Element::dimension(), 1>::Ones().eval();

    auto acceleration = Eigen::Matrix<double, Element::size(),
                                      Element::degrees_of_freedom()>::Zero()
                            .eval();
    for (std::size_t node = 0; node < Element::size(); node++)
      acceleration.template block<1, Element::dimension()>(node, 0) +=
          rigid_body_translational_acceleration;

    const auto mass_matrix = compute_mass_matrix(element);
    const auto force =
        (mass_matrix * acceleration.template reshaped<Eigen::RowMajor>(
                           Element::size() * Element::degrees_of_freedom(), 1))
            .template reshaped<Eigen::RowMajor>(Element::size(),
                                                Element::degrees_of_freedom());

    auto total_force =
        Eigen::Matrix<double, Element::dimension(), 1>::Zero().eval();
    for (std::size_t node = 0; node < Element::size(); node++)
      total_force += force.template block<1, Element::dimension()>(node, 0);

    EXPECT_THAT(
        total_force.isApprox(mass * rigid_body_translational_acceleration),
        true);
  }
};

TYPED_TEST_CASE_P(TimoshenkoBeamElementWithMass_Test);
TYPED_TEST_P(TimoshenkoBeamElementWithMass_Test,
             mass_matrix_is_semi_positive_definite) {
  this->check_semi_positive_definite();
}
TYPED_TEST_P(TimoshenkoBeamElementWithMass_Test,
             mass_matrix_complies_with_newtons_second_law) {
  this->check_second_law_of_motion();
}
REGISTER_TYPED_TEST_CASE_P(TimoshenkoBeamElementWithMass_Test,
                           mass_matrix_is_semi_positive_definite,
                           mass_matrix_complies_with_newtons_second_law);

template <std::size_t Dimension>
using TimoshenkoBeamElementWithMass =
    ElementWithMass<TimoshenkoBeamElement<Dimension>>;

using MassConfigurations =
    Types<ConfigurationWithMass<TimoshenkoBeamElementWithMass<2>,
                                &create_reference_element_axis<2>,
                                &timoshenko_beam_lumped_mass_matrix<2>>,
          ConfigurationWithMass<TimoshenkoBeamElementWithMass<3>,
                                &create_reference_element_axis<3>,
                                &timoshenko_beam_lumped_mass_matrix<3>>,
          ConfigurationWithMass<TimoshenkoBeamElementWithMass<2>,
                                &create_rotated_and_stretched_element_axis<2>,
                                &timoshenko_beam_lumped_mass_matrix<2>>,
          ConfigurationWithMass<TimoshenkoBeamElementWithMass<3>,
                                &create_rotated_and_stretched_element_axis<3>,
                                &timoshenko_beam_lumped_mass_matrix<3>>,

          ConfigurationWithMass<TimoshenkoBeamElementWithMass<2>,
                                &create_reference_element_axis<2>,
                                &timoshenko_beam_consistent_mass_matrix<2>>,
          ConfigurationWithMass<TimoshenkoBeamElementWithMass<3>,
                                &create_reference_element_axis<3>,
                                &timoshenko_beam_consistent_mass_matrix<3>>,
          ConfigurationWithMass<TimoshenkoBeamElementWithMass<2>,
                                &create_rotated_and_stretched_element_axis<2>,
                                &timoshenko_beam_consistent_mass_matrix<2>>,
          ConfigurationWithMass<TimoshenkoBeamElementWithMass<3>,
                                &create_rotated_and_stretched_element_axis<3>,
                                &timoshenko_beam_consistent_mass_matrix<3>>>;

INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElementWithMass_Test,
                              TimoshenkoBeamElementWithMass_Test,
                              MassConfigurations);
} // namespace
} // namespace elements
} // namespace ae108