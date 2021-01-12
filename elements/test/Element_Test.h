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

#pragma once

#include "ae108/elements/ElementBase.h"
#include "ae108/elements/automatic_forces.h"
#include "ae108/elements/automatic_stiffness_matrix.h"
#include "ae108/elements/compute_forces.h"
#include "ae108/elements/compute_stiffness_matrix.h"
#include "ae108/elements/tensor/as_matrix_of_columns.h"
#include <gmock/gmock.h>
#include <limits>
#include <type_traits>

namespace ae108 {
namespace elements {
namespace {

/**
 * @brief Compares the forces and the stiffness matrix to the numerical
 * approximation at example displacements.
 */
template <typename TestConfiguration> struct Element_Test : ::testing::Test {
  using Element = typename TestConfiguration::Element;
  using NodalDisplacements = typename Element::NodalDisplacements;

  static_assert(
      std::is_base_of<ElementBase<Element, typename Element::size_type,
                                  typename Element::value_type, Element::size(),
                                  Element::degrees_of_freedom()>,
                      Element>::value,
      "The element must derive from ElementBase.");

  Element element = TestConfiguration::create_element();
  const typename Element::Time time = TestConfiguration::create_time();

  /**
   * @brief Returns zero nodal displacements.
   */
  static NodalDisplacements zero_displacements() noexcept {
    return NodalDisplacements();
  }

  /**
   * @brief Returns nodal displacements equal to the identity matrix.
   */
  static NodalDisplacements identity_displacements() noexcept {
    NodalDisplacements displacements;
    auto matrix = tensor::as_matrix_of_columns(&displacements);
    matrix = matrix.Identity().eval();
    return displacements;
  }

  /**
   * @brief Returns nodal displacements where each entry is equal to 1/2.
   */
  static NodalDisplacements constant_displacements() noexcept {
    NodalDisplacements displacements;
    auto matrix = tensor::as_matrix_of_columns(&displacements);
    matrix = matrix.Constant(.5).eval();
    return displacements;
  }

  /**
   * @brief Returns nodal displacements where each upper diagonal entry is
   * equal to 1/2.
   */
  static NodalDisplacements upper_diagonal_displacements() noexcept {
    NodalDisplacements displacements = NodalDisplacements();
    auto matrix = tensor::as_matrix_of_columns(&displacements);
    matrix.template triangularView<Eigen::Upper>() = matrix.Constant(.5).eval();
    return displacements;
  }

  /**
   * @brief Returns nodal displacements where each lower diagonal entry is
   * equal to 1/2.
   */
  static NodalDisplacements lower_diagonal_displacements() noexcept {
    NodalDisplacements displacements = NodalDisplacements();
    auto matrix = tensor::as_matrix_of_columns(&displacements);
    matrix.template triangularView<Eigen::Lower>() = matrix.Constant(.5).eval();
    return displacements;
  }

  /**
   * @brief Returns random nodal displacements.
   */
  static NodalDisplacements random_displacements() noexcept {
    NodalDisplacements displacements;
    auto matrix = tensor::as_matrix_of_columns(&displacements);
    matrix = matrix.Random().eval();
    return displacements;
  }

  /**
   * @brief Computes the forces using the element.
   */
  typename Element::Forces
  computed_forces(const NodalDisplacements &displacements) const noexcept {
    return compute_forces(element, displacements, this->time);
  }

  /**
   * @brief Computes the forces using numerical differentiation.
   */
  typename Element::Forces
  approximated_forces(const NodalDisplacements &displacements) const noexcept {
    return automatic_forces(element, displacements, this->time);
  }

  /**
   * @brief Computes the tangent matrix using the element.
   */
  typename Element::StiffnessMatrix
  computed_stiffness_matrix(const NodalDisplacements &displacements) const
      noexcept {
    return compute_stiffness_matrix(element, displacements, this->time);
  }

  /**
   * @brief Computes the stress using numerical differentiation.
   */
  typename Element::StiffnessMatrix
  approximated_stiffness_matrix(const NodalDisplacements &displacements) const
      noexcept {
    return automatic_stiffness_matrix(element, displacements, this->time);
  }

  /**
   * @brief Returns true if and only if value is close to reference.
   */
  template <class ValueType, int Rows, int Cols>
  static bool is_approximately(
      const Eigen::Matrix<ValueType, Rows, Cols> &value,
      const Eigen::Matrix<ValueType, Rows, Cols> &reference) noexcept {
    constexpr auto epsilon =
        std::numeric_limits<typename Element::value_type>::epsilon();
    return (value.norm() <= epsilon && reference.norm() <= epsilon) ||
           value.isApprox(reference);
  }

  /**
   * Checks that the computed and the approximated forces are close.
   */
  void check_forces(const NodalDisplacements &displacements) const noexcept {
    const auto result = computed_forces(displacements);
    const auto approximation = approximated_forces(displacements);

    EXPECT_TRUE(
        is_approximately(tensor::as_matrix_of_columns(&result).eval(),
                         tensor::as_matrix_of_columns(&approximation).eval()));
  }

  /**
   * Checks that the computed and the approximated stiffness matrices are close.
   */
  void check_stiffness_matrix(const NodalDisplacements &displacements) const
      noexcept {
    const auto result = computed_stiffness_matrix(displacements);
    const auto approximation = approximated_stiffness_matrix(displacements);

    EXPECT_TRUE(is_approximately(result, approximation));
  }
};

TYPED_TEST_CASE_P(Element_Test);

TYPED_TEST_P(Element_Test,
             forces_are_derivative_of_energy_at_zero_displacements) {
  this->check_forces(this->zero_displacements());
}

TYPED_TEST_P(Element_Test,
             forces_are_derivative_of_energy_at_identity_displacements) {
  this->check_forces(this->identity_displacements());
}

TYPED_TEST_P(Element_Test,
             forces_are_derivative_of_energy_at_constant_displacements) {
  this->check_forces(this->constant_displacements());
}

TYPED_TEST_P(Element_Test,
             forces_are_derivative_of_energy_at_upper_diagonal_displacements) {
  this->check_forces(this->upper_diagonal_displacements());
}

TYPED_TEST_P(Element_Test,
             forces_are_derivative_of_energy_at_lower_diagonal_displacements) {
  this->check_forces(this->lower_diagonal_displacements());
}

TYPED_TEST_P(Element_Test,
             forces_are_derivative_of_energy_at_random_displacements) {
  this->check_forces(this->random_displacements());
}

TYPED_TEST_P(Element_Test,
             stiffness_is_derivative_of_stress_at_zero_displacements) {
  this->check_stiffness_matrix(this->zero_displacements());
}

TYPED_TEST_P(Element_Test,
             stiffness_is_derivative_of_stress_at_identity_displacements) {
  this->check_stiffness_matrix(this->identity_displacements());
}

TYPED_TEST_P(Element_Test,
             stiffness_is_derivative_of_stress_at_constant_displacements) {
  this->check_stiffness_matrix(this->constant_displacements());
}

TYPED_TEST_P(
    Element_Test,
    stiffness_is_derivative_of_stress_at_upper_diagonal_displacements) {
  this->check_stiffness_matrix(this->upper_diagonal_displacements());
}

TYPED_TEST_P(
    Element_Test,
    stiffness_is_derivative_of_stress_at_lower_diagonal_displacements) {
  this->check_stiffness_matrix(this->lower_diagonal_displacements());
}

TYPED_TEST_P(Element_Test,
             stiffness_is_derivative_of_stress_at_random_displacements) {
  this->check_stiffness_matrix(this->random_displacements());
}

REGISTER_TYPED_TEST_CASE_P(
    Element_Test, forces_are_derivative_of_energy_at_zero_displacements,
    forces_are_derivative_of_energy_at_identity_displacements,
    forces_are_derivative_of_energy_at_constant_displacements,
    forces_are_derivative_of_energy_at_upper_diagonal_displacements,
    forces_are_derivative_of_energy_at_lower_diagonal_displacements,
    forces_are_derivative_of_energy_at_random_displacements,
    stiffness_is_derivative_of_stress_at_zero_displacements,
    stiffness_is_derivative_of_stress_at_identity_displacements,
    stiffness_is_derivative_of_stress_at_constant_displacements,
    stiffness_is_derivative_of_stress_at_upper_diagonal_displacements,
    stiffness_is_derivative_of_stress_at_lower_diagonal_displacements,
    stiffness_is_derivative_of_stress_at_random_displacements);

} // namespace
} // namespace elements
} // namespace ae108