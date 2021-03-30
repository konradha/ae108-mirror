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
Eigen::Matrix<
    ValueType_, 2 * (Dimension_ + (Dimension_ * (Dimension_ - 1)) / 2),
    2 * (Dimension_ + (Dimension_ * (Dimension_ - 1)) / 2), Eigen::RowMajor>
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

  const double Phi_y = 12 * E * I_z * k_y / A / G / L / L;

  auto K = Eigen::Matrix<double, 12, 12, Eigen::RowMajor>::Zero().eval();

  const auto X = A * E / L;
  const auto Y1 = 12 * E * I_z / (1 + Phi_y) / L / L / L;
  const auto Y2 = 6 * E * I_z / (1 + Phi_y) / L / L;
  const auto Y3 = (4 + Phi_y) * E * I_z / (1 + Phi_y) / L;
  const auto Y4 = (2 - Phi_y) * E * I_z / (1 + Phi_y) / L;

  K(0, 0) = K(6, 6) = X;
  K(0, 6) = K(6, 0) = -X;
  K(1, 1) = K(7, 7) = Y1;
  K(1, 7) = K(7, 1) = -Y1;
  K(1, 5) = K(5, 1) = Y2;
  K(1, 11) = K(11, 1) = Y2;
  K(5, 7) = K(7, 5) = -Y2;
  K(7, 11) = K(11, 7) = -Y2;
  K(5, 5) = K(11, 11) = Y3;
  K(5, 11) = K(11, 5) = Y4;

  const auto I_y = properties.area_moment_y;
  const auto J_x = properties.polar_moment_x;
  const auto k_z = properties.shear_correction_factor_z;

  const double Phi_z = 12 * E * I_y * k_z / A / G / L / L;

  const auto Z1 = 12 * E * I_y / (1 + Phi_z) / L / L / L;
  const auto Z2 = 6 * E * I_y / (1 + Phi_z) / L / L;
  const auto Z3 = (4 + Phi_z) * E * I_y / (1 + Phi_z) / L;
  const auto Z4 = (2 - Phi_z) * E * I_y / (1 + Phi_z) / L;
  const auto S = G * J_x / L;

  K(2, 2) = K(8, 8) = Z1;
  K(2, 8) = K(8, 2) = -Z1;
  K(2, 4) = K(4, 2) = -Z2;
  K(2, 10) = K(10, 2) = -Z2;
  K(4, 8) = K(8, 4) = Z2;
  K(8, 10) = K(10, 8) = Z2;
  K(4, 4) = K(10, 10) = Z3;
  K(4, 10) = K(10, 4) = Z4;
  K(3, 3) = K(9, 9) = S;
  K(3, 9) = K(9, 3) = -S;

  return K * properties.weight;
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

  const double Phi_y = 12 * E * I_z * k_y / A / G / L / L;

  auto K = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Zero().eval();

  const auto X = A * E / L;
  const auto Y1 = 12 * E * I_z / (1 + Phi_y) / L / L / L;
  const auto Y2 = 6 * E * I_z / (1 + Phi_y) / L / L;
  const auto Y3 = (4 + Phi_y) * E * I_z / (1 + Phi_y) / L;
  const auto Y4 = (2 - Phi_y) * E * I_z / (1 + Phi_y) / L;

  K(0, 0) = K(3, 3) = X;
  K(0, 3) = K(3, 0) = -X;
  K(1, 1) = K(4, 4) = Y1;
  K(1, 4) = K(4, 1) = -Y1;
  K(1, 2) = K(2, 1) = Y2;
  K(1, 5) = K(5, 1) = Y2;
  K(2, 4) = K(4, 2) = -Y2;
  K(4, 5) = K(5, 4) = -Y2;
  K(2, 2) = K(5, 5) = Y3;
  K(2, 5) = K(5, 2) = Y4;

  return K * properties.weight;
}

template <class ValueType_, std::size_t Dimension_>
Eigen::Matrix<
    ValueType_, 2 * (Dimension_ + (Dimension_ * (Dimension_ - 1)) / 2),
    2 * (Dimension_ + (Dimension_ * (Dimension_ - 1)) / 2), Eigen::RowMajor>
rotation_matrix(const tensor::Tensor<ValueType_, Dimension_> &beam_orientation);

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.32
template <>
inline Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
rotation_matrix<double, 3>(const tensor::Tensor<double, 3> &beam_orientation) {

  const auto ax = tensor::as_vector(&beam_orientation) /
                  tensor::as_vector(&beam_orientation).norm();

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
rotation_matrix<double, 2>(const tensor::Tensor<double, 2> &beam_orientation) {

  const auto ax = tensor::as_vector(&beam_orientation) /
                  tensor::as_vector(&beam_orientation).norm();

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
                  Dimension_ + (Dimension_ * (Dimension_ - 1)) / 2> {
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