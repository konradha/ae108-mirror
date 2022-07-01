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

#include "ae108/elements/materialmodels/MaterialModelBase.h"
#include "ae108/elements/materialmodels/automatic_stress.h"
#include "ae108/elements/materialmodels/automatic_tangent_matrix.h"
#include "ae108/elements/materialmodels/compute_stress.h"
#include "ae108/elements/materialmodels/compute_tangent_matrix.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include "ae108/elements/tensor/as_two_tensor.h"
#include <gmock/gmock.h>
#include <limits>
#include <type_traits>

namespace ae108 {
namespace elements {
namespace materialmodels {
namespace {

/**
 * @brief Compares the stresses and the tangent matrix to the numerical
 * approximation at four example gradients.
 */
template <typename TestConfiguration>
struct MaterialModel_Test : ::testing::Test {
  using MaterialModel = typename TestConfiguration::MaterialModel;
  using DisplacementGradient = typename MaterialModel::DisplacementGradient;

  static_assert(
      std::is_base_of<MaterialModelBase<typename MaterialModel::size_type,
                                        typename MaterialModel::value_type,
                                        typename MaterialModel::real_type,
                                        MaterialModel::dimension()>,
                      MaterialModel>::value,
      "The material model must derive from MaterialModelBase.");

  MaterialModel model = TestConfiguration::create_model();
  const typename MaterialModel::Time time = TestConfiguration::create_time();

  /**
   * @brief Returns a displacement gradient of zeroes.
   */
  static DisplacementGradient zero_gradient() noexcept {
    return DisplacementGradient();
  }

  /**
   * @brief Returns a displacement gradient equal to the identity matrix.
   */
  static DisplacementGradient identity_gradient() noexcept {
    DisplacementGradient gradient;
    auto matrix = tensor::as_matrix_of_rows(&gradient);
    matrix = matrix.Identity().eval();
    return gradient;
  }

  /**
   * @brief Returns a displacement gradient where each entry is equal to 1/2.
   */
  static DisplacementGradient constant_gradient() noexcept {
    DisplacementGradient gradient;
    auto matrix = tensor::as_matrix_of_rows(&gradient);
    matrix = matrix.Constant(.5).eval();
    return gradient;
  }

  /**
   * @brief Returns a displacement gradient where each upper diagonal entry is
   * equal to 1/2.
   */
  static DisplacementGradient upper_diagonal_gradient() noexcept {
    DisplacementGradient gradient = DisplacementGradient();
    auto matrix = tensor::as_matrix_of_rows(&gradient);
    matrix.template triangularView<Eigen::Upper>() = matrix.Constant(.5).eval();
    return gradient;
  }

  /**
   * @brief Returns a displacement gradient where each lower diagonal entry is
   * equal to 1/2.
   */
  static DisplacementGradient lower_diagonal_gradient() noexcept {
    DisplacementGradient gradient = DisplacementGradient();
    auto matrix = tensor::as_matrix_of_rows(&gradient);
    matrix.template triangularView<Eigen::Lower>() = matrix.Constant(.5).eval();
    return gradient;
  }

  /**
   * @brief Returns a random displacement gradient.
   */
  static DisplacementGradient random_gradient() noexcept {
    DisplacementGradient gradient;
    auto matrix = tensor::as_matrix_of_rows(&gradient);
    matrix = matrix.Random().eval();
    return gradient;
  }

  /**
   * @brief Computes the stress using the material model.
   */
  typename MaterialModel::Stress
  computed_stress(const DisplacementGradient &gradient) const noexcept {
    return compute_stress(model, MaterialModel::unknown_id(), gradient,
                          this->time);
  }

  /**
   * @brief Computes the stress using numerical differentiation.
   */
  typename MaterialModel::Stress
  approximated_stress(const DisplacementGradient &gradient) const noexcept {
    return automatic_stress(model, MaterialModel::unknown_id(), gradient,
                            this->time);
  }

  /**
   * @brief Computes the tangent matrix using the material model.
   */
  typename MaterialModel::TangentMatrix
  computed_tangent_matrix(const DisplacementGradient &gradient) const noexcept {
    return compute_tangent_matrix(model, MaterialModel::unknown_id(), gradient,
                                  this->time);
  }

  /**
   * @brief Computes the stress using numerical differentiation.
   */
  typename MaterialModel::TangentMatrix approximated_tangent_matrix(
      const DisplacementGradient &gradient) const noexcept {
    return automatic_tangent_matrix(model, MaterialModel::unknown_id(),
                                    gradient, this->time);
  }

  /**
   * @brief Returns true if and only if value is close to reference.
   */
  template <class ValueType, int Rows, int Cols>
  static bool is_approximately(
      const Eigen::Matrix<ValueType, Rows, Cols> &value,
      const Eigen::Matrix<ValueType, Rows, Cols> &reference) noexcept {
    constexpr auto epsilon =
        std::numeric_limits<typename MaterialModel::real_type>::epsilon();
    return (value.norm() <= epsilon && reference.norm() <= epsilon) ||
           value.isApprox(reference);
  }

  /**
   * Checks that the computed and the approximated stresses are close.
   */
  void check_stresses(const DisplacementGradient &gradient) const noexcept {
    const auto result = computed_stress(gradient);
    const auto approximation = approximated_stress(gradient);

    EXPECT_TRUE(
        is_approximately(tensor::as_matrix_of_rows(&result).eval(),
                         tensor::as_matrix_of_rows(&approximation).eval()));
  }

  /**
   * Checks that the computed and the approximated tangent matrices are close.
   */
  void
  check_tangent_matrix(const DisplacementGradient &gradient) const noexcept {
    const auto result = computed_tangent_matrix(gradient);
    const auto approximation = approximated_tangent_matrix(gradient);

    EXPECT_TRUE(is_approximately(tensor::as_two_tensor(&result).eval(),
                                 tensor::as_two_tensor(&approximation).eval()));
  }
};

TYPED_TEST_CASE_P(MaterialModel_Test);

TYPED_TEST_P(MaterialModel_Test,
             stress_is_derivative_of_energy_at_zero_gradient) {
  this->check_stresses(this->zero_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             stress_is_derivative_of_energy_at_identity_gradient) {
  this->check_stresses(this->identity_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             stress_is_derivative_of_energy_at_constant_gradient) {
  this->check_stresses(this->constant_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             stress_is_derivative_of_energy_at_upper_diagonal_gradient) {
  this->check_stresses(this->upper_diagonal_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             stress_is_derivative_of_energy_at_lower_diagonal_gradient) {
  this->check_stresses(this->lower_diagonal_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             stress_is_derivative_of_energy_at_random_gradient) {
  this->check_stresses(this->random_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             tangent_is_derivative_of_stress_at_zero_gradient) {
  this->check_tangent_matrix(this->zero_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             tangent_is_derivative_of_stress_at_identity_gradient) {
  this->check_tangent_matrix(this->identity_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             tangent_is_derivative_of_stress_at_constant_gradient) {
  this->check_tangent_matrix(this->constant_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             tangent_is_derivative_of_stress_at_upper_diagonal_gradient) {
  this->check_tangent_matrix(this->upper_diagonal_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             tangent_is_derivative_of_stress_at_lower_diagonal_gradient) {
  this->check_tangent_matrix(this->lower_diagonal_gradient());
}

TYPED_TEST_P(MaterialModel_Test,
             tangent_is_derivative_of_stress_at_random_gradient) {
  this->check_tangent_matrix(this->random_gradient());
}

REGISTER_TYPED_TEST_CASE_P(
    MaterialModel_Test, stress_is_derivative_of_energy_at_zero_gradient,
    stress_is_derivative_of_energy_at_identity_gradient,
    stress_is_derivative_of_energy_at_constant_gradient,
    stress_is_derivative_of_energy_at_upper_diagonal_gradient,
    stress_is_derivative_of_energy_at_lower_diagonal_gradient,
    stress_is_derivative_of_energy_at_random_gradient,
    tangent_is_derivative_of_stress_at_zero_gradient,
    tangent_is_derivative_of_stress_at_identity_gradient,
    tangent_is_derivative_of_stress_at_constant_gradient,
    tangent_is_derivative_of_stress_at_upper_diagonal_gradient,
    tangent_is_derivative_of_stress_at_lower_diagonal_gradient,
    tangent_is_derivative_of_stress_at_random_gradient);

} // namespace
} // namespace materialmodels
} // namespace elements
} // namespace ae108