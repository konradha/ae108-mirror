// © 2021 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/elements/ElementWithMass.h"
#include "ae108/elements/TimoshenkoBeamElement.h"
#include "ae108/elements/compute_mass_matrix.h"
#include <Eigen/Dense>
#include <gmock/gmock.h>
#include <numeric>

using testing::DoubleEq;
using testing::Eq;
using testing::Not;
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

template <std::size_t Dimension>
using create_properties_t = decltype(&create_element_properties<Dimension>);

template <class Element_, create_axis_t<Element_::dimension()> CreateAxis_,
          create_properties_t<Element_::dimension()> CreateProperties_,
          compute_mass_t<Element_::dimension()> ComputeMass_>
struct ConfigurationWithMass {
  using Element = Element_;
  static Element create_element() noexcept {
    return Element(
        typename Element::Element(Element::StiffnessMatrix::Zero()),
        ComputeMass_(CreateAxis_(), CreateProperties_(),
                     create_element_density<Element_::dimension()>()));
  }
  static typename Element::Time calculate_mass() noexcept {
    const auto element_axis = CreateAxis_();
    const auto properties = create_element_properties<Element_::dimension()>();
    const auto density = create_element_density<Element_::dimension()>();
    return tensor::as_vector(&element_axis).norm() * properties.area * density;
  }
  static double calculate_length() noexcept {
    const auto axis = CreateAxis_();
    return tensor::as_vector(&axis).norm();
  }
  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

template <std::size_t Dimension>
using TimoshenkoBeamElementWithMass =
    ElementWithMass<TimoshenkoBeamElement<Dimension>>;

using ConfigurationsWithMass = Types<
    ConfigurationWithMass<
        TimoshenkoBeamElementWithMass<2>, &create_reference_element_axis<2>,
        &create_element_properties<2>, &timoshenko_beam_lumped_mass_matrix<2>>,
    ConfigurationWithMass<
        TimoshenkoBeamElementWithMass<3>, &create_reference_element_axis<3>,
        &create_element_properties<3>, &timoshenko_beam_lumped_mass_matrix<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElementWithMass_Test, Element_Test,
                              ConfigurationsWithMass);

MATCHER_P(IsEigenApprox, value,
          std::string(negation ? "not " : "") +
              "approximately equal to: " + ::testing::PrintToString(value)) {
  return arg.isApprox(value);
}

template <typename TestConfiguration>
struct TimoshenkoBeamElementWithMass_Test : ::testing::Test {
  using Element = typename TestConfiguration::Element;

  const Element element = TestConfiguration::create_element();
  const double mass = TestConfiguration::calculate_mass();
  const double length = TestConfiguration::calculate_length();

  /**
   * Checks that the mass matrix is positive (semi-)definite.
   */
  void check_semi_positive_definite() const noexcept {
    const auto mass_matrix = compute_mass_matrix(element);
    EXPECT_THAT(mass_matrix, IsEigenApprox(mass_matrix.transpose()));

    Eigen::LDLT<typename Element::MassMatrix> ldlt_of_mass_matrix(mass_matrix);
    EXPECT_THAT(ldlt_of_mass_matrix.info(), Not(Eq(Eigen::NumericalIssue)));
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

    EXPECT_THAT(total_force,
                IsEigenApprox(mass * rigid_body_translational_acceleration));
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

using MassConfigurations = Types<
    ConfigurationWithMass<
        TimoshenkoBeamElementWithMass<2>, &create_reference_element_axis<2>,
        &create_element_properties<2>, &timoshenko_beam_lumped_mass_matrix<2>>,
    ConfigurationWithMass<
        TimoshenkoBeamElementWithMass<3>, &create_reference_element_axis<3>,
        &create_element_properties<3>, &timoshenko_beam_lumped_mass_matrix<3>>,
    ConfigurationWithMass<TimoshenkoBeamElementWithMass<2>,
                          &create_rotated_and_stretched_element_axis<2>,
                          &create_element_properties<2>,
                          &timoshenko_beam_lumped_mass_matrix<2>>,
    ConfigurationWithMass<TimoshenkoBeamElementWithMass<3>,
                          &create_rotated_and_stretched_element_axis<3>,
                          &create_element_properties<3>,
                          &timoshenko_beam_lumped_mass_matrix<3>>,

    ConfigurationWithMass<TimoshenkoBeamElementWithMass<2>,
                          &create_reference_element_axis<2>,
                          &create_element_properties<2>,
                          &timoshenko_beam_consistent_mass_matrix<2>>,
    ConfigurationWithMass<TimoshenkoBeamElementWithMass<3>,
                          &create_reference_element_axis<3>,
                          &create_element_properties<3>,
                          &timoshenko_beam_consistent_mass_matrix<3>>,
    ConfigurationWithMass<TimoshenkoBeamElementWithMass<2>,
                          &create_rotated_and_stretched_element_axis<2>,
                          &create_element_properties<2>,
                          &timoshenko_beam_consistent_mass_matrix<2>>,
    ConfigurationWithMass<TimoshenkoBeamElementWithMass<3>,
                          &create_rotated_and_stretched_element_axis<3>,
                          &create_element_properties<3>,
                          &timoshenko_beam_consistent_mass_matrix<3>>>;

INSTANTIATE_TYPED_TEST_CASE_P(TimoshenkoBeamElementWithMass_Test,
                              TimoshenkoBeamElementWithMass_Test,
                              MassConfigurations);

template <std::size_t Dimension_>
TimoshenkoBeamProperties<double, Dimension_>
create_euler_bernoulli_properties() noexcept {
  static_assert(Dimension_ == 2 || Dimension_ == 3,
                "Only dimensions 2 and 3 are supported.");
  if constexpr (Dimension_ == 2) {
    return {1., 1., 1., 1., 0.};
  } else {
    return {1., 1., 1., 1., 1., 0., 0., 0.};
  }
}

struct TimoshenkoBeamElement2D_EulerBernoulliTest
    : TimoshenkoBeamElementWithMass_Test<ConfigurationWithMass<
          TimoshenkoBeamElementWithMass<2>, &create_reference_element_axis<2>,
          &create_euler_bernoulli_properties<2>,
          &timoshenko_beam_consistent_mass_matrix<2>>> {};

TEST_F(TimoshenkoBeamElement2D_EulerBernoulliTest,
       special_case_is_equal_to_euler_bernoulli_mass_matrix) {
  const auto matrix = compute_mass_matrix(element);

  // see Felippa et al (2015), "Mass Matrix Templates: General Description
  // and 1D Examples", eq. 148, http://dx.doi.org/10.1007/s11831-014-9108-x
  //
  // Note that degrees of freedom 0, 3 represent the displacement in axial
  // direction.
  const auto reference = [&]() {
    auto reference = decltype(matrix)::Zero().eval();
    reference(0, 0) = 420. / 3.;
    reference(0, 3) = reference(3, 0) = 420 / 6.;
    reference(3, 3) = 420. / 3.;
    reference(1, 1) = 156.;
    reference(1, 2) = reference(2, 1) = 22. * length;
    reference(1, 4) = reference(4, 1) = 54;
    reference(1, 5) = reference(5, 1) = -13. * length;
    reference(2, 2) = 4. * length * length;
    reference(2, 4) = reference(4, 2) = 13. * length;
    reference(2, 5) = reference(5, 2) = -3. * length * length;
    reference(4, 4) = 156.;
    reference(4, 5) = reference(5, 4) = -22. * length;
    reference(5, 5) = 4 * length * length;
    return ((mass / 420.) * reference).eval();
  }();

  EXPECT_THAT(matrix, IsEigenApprox(reference));
}

struct TimoshenkoBeamElement3D_EulerBernoulliTest
    : TimoshenkoBeamElementWithMass_Test<ConfigurationWithMass<
          TimoshenkoBeamElementWithMass<3>, &create_reference_element_axis<3>,
          &create_euler_bernoulli_properties<3>,
          &timoshenko_beam_consistent_mass_matrix<3>>> {};

TEST_F(TimoshenkoBeamElement3D_EulerBernoulliTest,
       special_case_is_equal_to_euler_bernoulli_mass_matrix) {
  const auto matrix = compute_mass_matrix(element);

  // see Felippa et al (2015), "Mass Matrix Templates: General Description
  // and 1D Examples", eq. 148, http://dx.doi.org/10.1007/s11831-014-9108-x
  //
  // Note that degrees of freedom 0, 6 represent the displacement in axial
  // direction. Degrees of freedom 3, 9 represent torsion around the beam axis.
  const auto reference = [&]() {
    auto reference = decltype(matrix)::Zero().eval();
    reference(0, 0) = reference(6, 6) = 420. / 3.;
    reference(0, 6) = reference(6, 0) = 420 / 6.;
    reference(1, 1) = reference(2, 2) = 156.;
    reference(1, 5) = reference(5, 1) = 22. * length;
    reference(2, 4) = reference(4, 2) = -22. * length;
    reference(1, 7) = reference(7, 1) = 54;
    reference(2, 8) = reference(8, 2) = 54;
    reference(1, 11) = reference(11, 1) = -13. * length;
    reference(2, 10) = reference(10, 2) = 13. * length;
    reference(4, 4) = reference(5, 5) = 4. * length * length;
    reference(4, 10) = reference(10, 4) = -3. * length * length;
    reference(5, 11) = reference(11, 5) = -3. * length * length;
    reference(7, 7) = reference(8, 8) = 156.;
    reference(7, 11) = reference(11, 7) = -22. * length;
    reference(8, 10) = reference(10, 8) = 22. * length;
    reference(4, 8) = reference(8, 4) = -13. * length;
    reference(7, 5) = reference(5, 7) = 13. * length;
    reference(11, 11) = reference(10, 10) = 4 * length * length;

    return ((mass / 420.) * reference).eval();
  }();

  EXPECT_THAT(matrix, IsEigenApprox(reference));
}

struct TimoshenkoBeamElement2D_LumpedMassTest
    : TimoshenkoBeamElementWithMass_Test<ConfigurationWithMass<
          TimoshenkoBeamElementWithMass<2>, &create_reference_element_axis<2>,
          &create_euler_bernoulli_properties<2>,
          &timoshenko_beam_lumped_mass_matrix<2>>> {};

TEST_F(TimoshenkoBeamElement2D_LumpedMassTest, lumped_mass_matrix_is_correct) {
  const auto matrix = compute_mass_matrix(element);

  // see Felippa et al (2015), "Mass Matrix Templates: General Description
  // and 1D Examples", eq. 123, http://dx.doi.org/10.1007/s11831-014-9108-x
  //
  // Note that degrees of freedom 0, 3 represent the displacement in axial
  // direction.
  const auto reference = [&]() {
    auto reference = decltype(matrix)::Zero().eval();
    reference(0, 0) = .5 * mass;
    reference(1, 1) = .5 * mass;
    reference(2, 2) = 1. / 24. * length * length * mass;
    reference(3, 3) = .5 * mass;
    reference(4, 4) = .5 * mass;
    reference(5, 5) = 1. / 24. * length * length * mass;
    return reference;
  }();

  EXPECT_THAT(matrix, IsEigenApprox(reference));
}

struct TimoshenkoBeamElement3D_LumpedMassTest
    : TimoshenkoBeamElementWithMass_Test<ConfigurationWithMass<
          TimoshenkoBeamElementWithMass<3>, &create_reference_element_axis<3>,
          &create_euler_bernoulli_properties<3>,
          &timoshenko_beam_lumped_mass_matrix<3>>> {};

TEST_F(TimoshenkoBeamElement3D_LumpedMassTest, lumped_mass_matrix_is_correct) {
  const auto matrix = compute_mass_matrix(element);

  // see Felippa et al (2015), "Mass Matrix Templates: General Description
  // and 1D Examples", eq. 123, http://dx.doi.org/10.1007/s11831-014-9108-x
  //
  // Note that degrees of freedom 0, 3 represent the displacement in axial
  // direction.
  const auto reference = [&]() {
    auto reference = decltype(matrix)::Zero().eval();
    reference(0, 0) = .5 * mass;
    reference(1, 1) = .5 * mass;
    reference(2, 2) = .5 * mass;
    reference(4, 4) = 1. / 24. * length * length * mass;
    reference(5, 5) = 1. / 24. * length * length * mass;
    reference(6, 6) = .5 * mass;
    reference(7, 7) = .5 * mass;
    reference(8, 8) = .5 * mass;
    reference(10, 10) = 1. / 24. * length * length * mass;
    reference(11, 11) = 1. / 24. * length * length * mass;
    return reference;
  }();

  EXPECT_THAT(matrix, IsEigenApprox(reference));
}

} // namespace
} // namespace elements
} // namespace ae108