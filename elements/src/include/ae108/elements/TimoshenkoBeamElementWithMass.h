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

template <class RealType_, std::size_t Dimension_>
struct TimoshenkoBeamWithMassProperties
    : TimoshenkoBeamProperties<RealType_, Dimension_> {
  RealType_ density;
};

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
    const TimoshenkoBeamWithMassProperties<double, Dimension_>
        &properties) noexcept;

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
    const TimoshenkoBeamWithMassProperties<double, Dimension_>
        &properties) noexcept;

template <class Element_>
struct ElementWithMass final
    : ElementBase<ElementWithMass<Element_>, typename Element_::size_type,
                  typename Element_::value_type, typename Element_::real_type,
                  Element_::size(), Element_::dimension(),
                  Element_::degrees_of_freedom()> {
public:
  using Element = Element_;
  using MassMatrix = typename ElementWithMass::StiffnessMatrix;

  explicit ElementWithMass(Element element, MassMatrix mass_matrix) noexcept
      : element_(std::move(element)), mass_matrix_(std::move(mass_matrix)) {}

  const MassMatrix &mass_matrix() const & { return mass_matrix_; }

private:
  Element element_;
  MassMatrix mass_matrix_;
};

template <class Element_>
struct ComputeEnergyTrait<ElementWithMass<Element_>>
    : ComputeEnergyTrait<Element_> {};

template <class Element_>
struct ComputeForcesTrait<ElementWithMass<Element_>>
    : ComputeForcesTrait<Element_> {};

template <class Element_>
struct ComputeStiffnessMatrixTrait<ElementWithMass<Element_>>
    : ComputeStiffnessMatrixTrait<Element_> {};

template <class Element_>
struct ComputeMassMatrixTrait<ElementWithMass<Element_>> {
  template <class Element>
  auto &operator()(const Element &element) const &noexcept {
    return element.mass_matrix();
  }
};
} // namespace elements
} // namespace ae108