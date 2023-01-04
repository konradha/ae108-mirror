// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#pragma once

#include "ae108/elements/materialmodels/AutomaticStressTrait.h"
#include "ae108/elements/materialmodels/AutomaticTangentMatrixTrait.h"
#include "ae108/elements/materialmodels/ComputeEnergyTrait.h"
#include "ae108/elements/materialmodels/ComputeStressTrait.h"
#include "ae108/elements/materialmodels/ComputeTangentMatrixTrait.h"
#include "ae108/elements/materialmodels/MaterialModelBase.h"

namespace ae108 {
namespace elements {
namespace materialmodels {

/**
 * @brief A minimal example of a material model. Note that only the energy is
 * defined.
 */
template <std::size_t Dimension_>
struct Minimal final
    : MaterialModelBase<std::size_t, double, double, Dimension_> {};

/**
 * @brief Always returns 0.
 */
template <std::size_t Dimension_>
struct ComputeEnergyTrait<Minimal<Dimension_>> {
  template <class MaterialModel>
  typename MaterialModel::Energy
  operator()(const MaterialModel &, const typename MaterialModel::size_type,
             const typename MaterialModel::DisplacementGradient &,
             const typename MaterialModel::Time) noexcept {
    return 0.;
  }
};

/**
 * @brief Computes the stress by differentiating the energy.
 */
template <std::size_t Dimension_>
struct ComputeStressTrait<Minimal<Dimension_>>
    : AutomaticStressTrait<Minimal<Dimension_>> {};

/**
 * @brief Computes the tangent matrix by differentiating the stress.
 */
template <std::size_t Dimension_>
struct ComputeTangentMatrixTrait<Minimal<Dimension_>>
    : AutomaticTangentMatrixTrait<Minimal<Dimension_>> {};

} // namespace materialmodels
} // namespace elements
} // namespace ae108
