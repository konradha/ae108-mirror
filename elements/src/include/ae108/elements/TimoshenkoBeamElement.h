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

  value_type density;

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

  value_type density;

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

  const tensor::Tensor<double, 12, 12> matrix = {{
      {{X, 0., 0., 0., 0., 0., -X, 0., 0., 0., 0., 0.}},
      {{_, Y1, 0., 0., 0., Y2, 0., -Y1, 0., 0., 0., Y2}},
      {{_, _, Z1, 0., -Z2, 0., 0., 0., -Z1, 0., -Z2, 0.}},
      {{_, _, _, S, 0., 0., 0., 0., 0., -S, 0., 0.}},
      {{_, _, _, _, Z3, 0., 0., 0., Z2, 0., Z4, 0.}},
      {{_, _, _, _, _, Y3, 0. - Y2, 0., 0., 0., Y4}},
      {{_, _, _, _, _, _, X, 0., 0., 0., 0., 0.}},
      {{_, _, _, _, _, _, _, Y1, 0., 0., 0., -Y2}},
      {{_, _, _, _, _, _, _, _, Z1, 0., Z2, 0.}},
      {{_, _, _, _, _, _, _, _, _, S, 0., 0.}},
      {{_, _, _, _, _, _, _, _, _, _, Z3, 0.}},
      {{_, _, _, _, _, _, _, _, _, _, _, Y3}},
  }};

  return properties.weight * tensor::as_matrix_of_rows(&matrix)
                                 .template selfadjointView<Eigen::Upper>();
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

  const tensor::Tensor<double, 6, 6> matrix = {{
      {{X, 0., 0., -X, 0., 0.}},
      {{0., Y1, Y2, 0., -Y1, Y2}},
      {{0., Y2, Y3, 0., -Y2, Y4}},
      {{-X, 0., 0., X, 0., 0.}},
      {{0, -Y1, -Y2, 0., Y1, -Y2}},
      {{0., Y2, Y4, 0., -Y2, Y3}},
  }};

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

  const auto ax =
      tensor::as_vector(&orientation) / tensor::as_vector(&orientation).norm();

  auto Lambda = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>::Zero().eval();
  Lambda.col(0) << ax(0), -ax(0) * ax(1), -ax(2);              // l1,l2,l3
  Lambda.col(1) << ax(1), (ax(0) * ax(0) + ax(2) * ax(2)), 0.; // m1,m2,m3
  Lambda.col(2) << ax(2), -ax(1) * ax(2), ax(0);               // n1,n2,n3
  Lambda.block(1, 0, 2, 3) /= std::sqrt(ax(0) * ax(0) + ax(2) * ax(2));

  if (fabs(ax(0)) < 1e-4 && fabs(ax(2)) < 1e-4)
    Lambda << 0, ax(1), 0, -ax(1), 0, 0, 0, 0, 1;

  auto T = Eigen::Matrix<double, 12, 12, Eigen::RowMajor>::Zero().eval();

  T.block(0, 0, 3, 3) = Lambda;
  T.block(3, 3, 3, 3) = Lambda;
  T.block(6, 6, 3, 3) = Lambda;
  T.block(9, 9, 3, 3) = Lambda;

  return T;
}

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.31
template <>
inline Eigen::Matrix<double, 6, 6, Eigen::RowMajor>
rotation_matrix<double, 2>(const tensor::Tensor<double, 2> &orientation) {

  const auto ax =
      tensor::as_vector(&orientation) / tensor::as_vector(&orientation).norm();

  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Lambda =
      Eigen::Matrix<double, 3, 3, Eigen::RowMajor>::Zero();

  Lambda.col(0) << ax(0), -ax(1), 0.;
  Lambda.col(1) << ax(1), ax(0), 0.;
  Lambda.col(2) << 0., 0., 1.;

  auto T = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Zero().eval();

  T.block(0, 0, 3, 3) = Lambda;
  T.block(3, 3, 3, 3) = Lambda;

  return T;
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