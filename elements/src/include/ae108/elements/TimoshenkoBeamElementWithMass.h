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

#include "ae108/elements/ComputeMassMatrixTrait.h"
#include "ae108/elements/TimoshenkoBeamElement.h"

namespace ae108 {
namespace elements {

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
    const TimoshenkoBeamProperties<double, Dimension_> &properties) noexcept;

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
    const TimoshenkoBeamProperties<double, Dimension_> &properties) noexcept;

template <std::size_t Dimension_, class ValueType_ = double,
          class RealType_ = double>
struct TimoshenkoBeamElementWithMass final
    : ElementBase<
          TimoshenkoBeamElementWithMass<Dimension_, ValueType_, RealType_>,
                  std::size_t, ValueType_, RealType_, 2,
                  (Dimension_ * (Dimension_ + 1)) / 2> {
public:
  using Element = TimoshenkoBeamElement<Dimension_>;
  using StiffnessMatrix =
      typename TimoshenkoBeamElementWithMass::StiffnessMatrix;
  using MassMatrix = StiffnessMatrix;

  explicit TimoshenkoBeamElementWithMass(Element element,
                                         MassMatrix mass_matrix) noexcept
      : element_(std::move(element)), mass_matrix_(std::move(mass_matrix)) {}

  const StiffnessMatrix &stiffness_matrix() const {
    return element_.stiffness_matrix();
  }

  static constexpr typename Element::size_type dimension() {
    return Dimension_;
  }

  const MassMatrix &mass_matrix() const { return mass_matrix_; }

  const MassMatrix computeMassMatrix() const { return mass_matrix_; }

private:
  TimoshenkoBeamElement<Dimension_> element_;
  MassMatrix mass_matrix_;
};

template <std::size_t Dimension_>
struct ComputeEnergyTrait<TimoshenkoBeamElementWithMass<Dimension_>>
    : ComputeEnergyTrait<TimoshenkoBeamElement<Dimension_>> {};

template <std::size_t Dimension_>
struct ComputeForcesTrait<TimoshenkoBeamElementWithMass<Dimension_>>
    : ComputeForcesTrait<TimoshenkoBeamElement<Dimension_>> {};

template <std::size_t Dimension_>
struct ComputeStiffnessMatrixTrait<TimoshenkoBeamElementWithMass<Dimension_>>
    : ComputeStiffnessMatrixTrait<TimoshenkoBeamElement<Dimension_>> {};

template <std::size_t Dimension_>
struct ComputeMassMatrixTrait<TimoshenkoBeamElementWithMass<Dimension_>> {
  typename TimoshenkoBeamElementWithMass<Dimension_>::MassMatrix operator()(
      const TimoshenkoBeamElementWithMass<Dimension_> &element) const noexcept {
    return element.mass_matrix();
  }
};
} // namespace elements
} // namespace ae108