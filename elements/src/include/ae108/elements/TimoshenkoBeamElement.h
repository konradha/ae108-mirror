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

template <class RealType_, std::size_t Dimension_>
struct TimoshenkoBeamProperties;

template <class RealType_> struct TimoshenkoBeamProperties<RealType_, 3> {
  using real_type = RealType_;

  real_type young_modulus;
  real_type shear_modulus;

  real_type shear_correction_factor_y;
  real_type shear_correction_factor_z;

  real_type area;
  real_type area_moment_z;
  real_type area_moment_y;
  real_type polar_moment_x;
};

template <class RealType_> struct TimoshenkoBeamProperties<RealType_, 2> {
  using real_type = RealType_;

  real_type young_modulus;
  real_type shear_modulus;

  real_type shear_correction_factor_y;

  real_type area;

  real_type area_moment_z;
};

/**
 * @brief Computes the stiffness matrix for a Timoshenko beam with the given
 * axis and the given properties.
 * @tparam Dimension_ The dimenion of the physical space. Only dimensions 2 and
 * 3 are supported.
 */
template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
timoshenko_beam_stiffness_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TimoshenkoBeamProperties<double, Dimension_> &properties) noexcept;

/**
 * @brief Computes the lumped mass matrix for a Timoshenko beam with the given
 * axis and the given properties.
 * @tparam Dimension_ The dimenion of the physical space. Only dimensions 2 and
 * 3 are supported.
 */
template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
timoshenko_beam_lumped_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TimoshenkoBeamProperties<double, Dimension_> &properties,
    const double density) noexcept;

/**
 * @brief Computes the consistent mass matrix for a Timoshenko beam with the
 * given axis and the given properties.
 * @tparam Dimension_ The dimenion of the physical space. Only dimensions 2 and
 * 3 are supported.
 */
template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
timoshenko_beam_consistent_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TimoshenkoBeamProperties<double, Dimension_> &properties,
    const double density) noexcept;

/**
 * @brief Implementation of the closed-form Timoshenko beam element as presented
 * in Cook et. al (2002), "Concepts and applications of Finite Element
 * Analysis", 4th ed., pp.24-32
 */
template <std::size_t Dimension_, class ValueType_ = double,
          class RealType_ = double>
struct TimoshenkoBeamElement final
    : ElementBase<TimoshenkoBeamElement<Dimension_, ValueType_, RealType_>,
                  std::size_t, ValueType_, RealType_, 2, Dimension_,
                  (Dimension_ * (Dimension_ + 1)) / 2> {
public:
  explicit TimoshenkoBeamElement(
      typename TimoshenkoBeamElement::StiffnessMatrix matrix) noexcept
      : stiffness_matrix_(std::move(matrix)) {}

  const typename TimoshenkoBeamElement::StiffnessMatrix &
  stiffness_matrix() const {
    return stiffness_matrix_;
  }

private:
  typename TimoshenkoBeamElement::StiffnessMatrix stiffness_matrix_;
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeEnergyTrait<
    TimoshenkoBeamElement<Dimension_, ValueType_, RealType_>> {
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

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeForcesTrait<
    TimoshenkoBeamElement<Dimension_, ValueType_, RealType_>> {
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

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeStiffnessMatrixTrait<
    TimoshenkoBeamElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &element,
             const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return element.stiffness_matrix();
  }
};

} // namespace elements
} // namespace ae108