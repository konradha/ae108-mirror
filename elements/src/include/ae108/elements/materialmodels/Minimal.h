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
struct Minimal final : MaterialModelBase<std::size_t, double, Dimension_> {};

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
