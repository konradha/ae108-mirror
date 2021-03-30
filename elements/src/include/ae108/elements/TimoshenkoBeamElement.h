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

#pragma once

#include "ae108/elements/ComputeEnergyTrait.h"
#include "ae108/elements/ComputeForcesTrait.h"
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include "ae108/elements/tensor/as_vector.h"
#include <Eigen/Geometry>

namespace ae108 {
namespace elements {
namespace timoshenko {

template <class ValueType_, std::size_t Dimension_> struct Properties;

template <class ValueType_> struct Properties<ValueType_, 3> {
  using value_type = ValueType_;

  value_type young_modulus;
  value_type shear_modulus;

  value_type shear_correction_factor_y;
  value_type shear_correction_factor_z;

  value_type area;

  value_type area_moment_z;
  value_type area_moment_y;
  value_type polar_moment_x;

  value_type weight;
};

template <class ValueType_> struct Properties<ValueType_, 2> {
  using value_type = ValueType_;

  value_type young_modulus;
  value_type shear_modulus;

  value_type shear_correction_factor_y;

  value_type area;

  value_type area_moment_z;

  value_type weight;
};

/**
 * @brief Computes the stiffness matrix of a reference beam with the given
 * properties.
 */
template <class ValueType_, std::size_t Dimension_>
Eigen::Matrix<ValueType_, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
stiffness_matrix(const Properties<ValueType_, Dimension_> &properties,
                 const ValueType_ length);

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.27
template <>
inline Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
stiffness_matrix<double, 3>(const Properties<double, 3> &properties,
                            const double length) {
  const auto L = length;
  const auto A = properties.area;
  const auto E = properties.young_modulus;
  const auto G = properties.shear_modulus;
  const auto I_z = properties.area_moment_z;
  const auto k_y = properties.shear_correction_factor_y;
  const auto I_y = properties.area_moment_y;
  const auto J_x = properties.polar_moment_x;
  const auto k_z = properties.shear_correction_factor_z;

  const double phi_y = 12 * E * I_z * k_y / A / G / L / L;
  const double phi_z = 12 * E * I_y * k_z / A / G / L / L;

  const auto X = A * E / L;
  const auto Y1 = 12 * E * I_z / (1 + phi_y) / L / L / L;
  const auto Y2 = 6 * E * I_z / (1 + phi_y) / L / L;
  const auto Y3 = (4 + phi_y) * E * I_z / (1 + phi_y) / L;
  const auto Y4 = (2 - phi_y) * E * I_z / (1 + phi_y) / L;
  const auto Z1 = 12 * E * I_y / (1 + phi_z) / L / L / L;
  const auto Z2 = 6 * E * I_y / (1 + phi_z) / L / L;
  const auto Z3 = (4 + phi_z) * E * I_y / (1 + phi_z) / L;
  const auto Z4 = (2 - phi_z) * E * I_y / (1 + phi_z) / L;
  const auto S = G * J_x / L;

  const auto _ = 0.;

  // clang-format off
  const tensor::Tensor<double, 12, 12> matrix = {{
      {{  X,   _,   _,   _,   _,   _,  -X,   _,   _,   _,   _,   _}},
      {{  _,  Y1,   _,   _,   _,  Y2,   _, -Y1,   _,   _,   _,  Y2}},
      {{  _,   _,  Z1,   _, -Z2,   _,   _,   _, -Z1,   _, -Z2,   _}},
      {{  _,   _,   _,   S,   _,   _,   _,   _,   _,  -S,   _,   _}},
      {{  _,   _, -Z2,   _,  Z3,   _,   _,   _,  Z2,   _,  Z4,   _}},
      {{  _,  Y2,   _,   _,   _,  Y3,   _, -Y2,   _,   _,   _,  Y4}},
      {{ -X,   _,   _,   _,   _,   _,   X,   _,   _,   _,   _,   _}},
      {{  _, -Y1,   _,   _,   _, -Y2,   _,  Y1,   _,   _,   _, -Y2}},
      {{  _,   _, -Z1,   _,  Z2,   _,   _,   _,  Z1,   _,  Z2,   _}},
      {{  _,   _,   _,  -S,   _,   _,   _,   _,   _,   S,   _,   _}},
      {{  _,   _, -Z2,   _,  Z4,   _,   _,   _,  Z2,   _,  Z3,   _}},
      {{  _,  Y2,   _,   _,   _,  Y4,   _, -Y2,   _,   _,   _,  Y3}},
  }};
  // clang-format on

  return properties.weight * tensor::as_matrix_of_rows(&matrix);
}

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.26
template <>
inline Eigen::Matrix<double, 6, 6, Eigen::RowMajor>
stiffness_matrix<double, 2>(const Properties<double, 2> &properties,
                            const double length) {
  const auto L = length;
  const auto A = properties.area;
  const auto E = properties.young_modulus;
  const auto G = properties.shear_modulus;
  const auto I_z = properties.area_moment_z;
  const auto k_y = properties.shear_correction_factor_y;

  const double phi_y = 12 * E * I_z * k_y / A / G / L / L;

  const auto X = A * E / L;
  const auto Y1 = 12 * E * I_z / (1 + phi_y) / L / L / L;
  const auto Y2 = 6 * E * I_z / (1 + phi_y) / L / L;
  const auto Y3 = (4 + phi_y) * E * I_z / (1 + phi_y) / L;
  const auto Y4 = (2 - phi_y) * E * I_z / (1 + phi_y) / L;

  const auto _ = 0.;

  // clang-format off
  const tensor::Tensor<double, 6, 6> matrix = {{
      {{  X,   _,   _,  -X,   _,   _}},
      {{  _,  Y1,  Y2,   _, -Y1,  Y2}},
      {{  _,  Y2,  Y3,   _, -Y2,  Y4}},
      {{ -X,   _,   _,   X,   _,   _}},
      {{  _, -Y1, -Y2,   _,  Y1, -Y2}},
      {{  _,  Y2,  Y4,   _, -Y2,  Y3}},
  }};
  // clang-format on

  return properties.weight * tensor::as_matrix_of_rows(&matrix);
}

template <class ValueType_, std::size_t Dimension_>
Eigen::Matrix<ValueType_, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
rotation_matrix(const tensor::Tensor<ValueType_, Dimension_> &orientation);

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.32
template <>
inline Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
rotation_matrix<double, 3>(const tensor::Tensor<double, 3> &orientation) {
  // rotation that maps the normalized orientation vector to (1, 0, 0)
  const auto Lambda = Eigen::Quaternion<double>()
                          .FromTwoVectors(tensor::as_vector(&orientation),
                                          Eigen::Vector3d::UnitX())
                          .normalized()
                          .toRotationMatrix();

  auto result = Eigen::Matrix<double, 12, 12, Eigen::RowMajor>::Zero().eval();

  result.block(0, 0, 3, 3) = Lambda;
  result.block(3, 3, 3, 3) = Lambda;
  result.block(6, 6, 3, 3) = Lambda;
  result.block(9, 9, 3, 3) = Lambda;

  return result;
}

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.31
template <>
inline Eigen::Matrix<double, 6, 6, Eigen::RowMajor>
rotation_matrix<double, 2>(const tensor::Tensor<double, 2> &orientation) {
  const auto normalized =
      tensor::as_vector(&orientation) / tensor::as_vector(&orientation).norm();

  // rotation that maps the normalized orientation vector to (1, 0)
  const tensor::Tensor<double, 2, 2> Lambda = {{
      {{normalized[0], normalized[1]}},
      {{-normalized[1], normalized[0]}},
  }};

  auto result = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Zero().eval();

  result.block(0, 0, 2, 2) = tensor::as_matrix_of_rows(&Lambda);
  result(2, 2) = 1.;
  result.block(3, 3, 2, 2) = tensor::as_matrix_of_rows(&Lambda);
  result(5, 5) = 1.;

  return result;
}

/**
 * @brief Implementation of the closed-form Timoshenko beam element as presented
 * in Cook et. al (2002), "Concepts and applications of Finite Element
 * Analysis", 4th ed., pp.24-32
 */
template <std::size_t Dimension_>
struct BeamElement final
    : ElementBase<BeamElement<Dimension_>, std::size_t, double, 2,
                  (Dimension_ * (Dimension_ + 1)) / 2> {
public:
  using value_type = typename BeamElement::value_type;
  using size_type = typename BeamElement::size_type;

  using Vector = tensor::Tensor<value_type, Dimension_>;

  explicit BeamElement(
      const Vector &element_axis,
      const Properties<value_type, Dimension_> &properties) noexcept {

    const auto reference_stiffness_matrix =
        timoshenko::stiffness_matrix<value_type, BeamElement::dimension()>(
            properties, tensor::as_vector(&element_axis).norm());

    const auto rotation_matrix =
        timoshenko::rotation_matrix<value_type, BeamElement::dimension()>(
            element_axis);

    stiffness_matrix_ = rotation_matrix.transpose() *
                        reference_stiffness_matrix * rotation_matrix;
  }

  const typename BeamElement::StiffnessMatrix &stiffness_matrix() const {
    return stiffness_matrix_;
  }

  static constexpr size_type dimension() { return Dimension_; }

private:
  typename BeamElement::StiffnessMatrix stiffness_matrix_;
};
} // namespace timoshenko

template <std::size_t Dimension_>
struct ComputeEnergyTrait<timoshenko::BeamElement<Dimension_>> {
  template <class Element>
  typename Element::Energy
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    const auto v = tensor::as_vector(&u);
    return typename Element::Energy{.5} * v.transpose() *
           element.stiffness_matrix() * v;
  }
};

template <std::size_t Dimension_>
struct ComputeForcesTrait<timoshenko::BeamElement<Dimension_>> {
  template <class Element>
  typename Element::Forces
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    typename Element::Forces forces;
    tensor::as_vector(&forces) =
        element.stiffness_matrix() * tensor::as_vector(&u);
    return forces;
  }
};

template <std::size_t Dimension_>
struct ComputeStiffnessMatrixTrait<timoshenko::BeamElement<Dimension_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &element,
             const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return element.stiffness_matrix();
  };
};

} // namespace elements
} // namespace ae108