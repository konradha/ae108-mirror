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

#include "ae108/elements/materialmodels/compute_energy.h"
#include "ae108/elements/tensor/differentiate.h"

namespace ae108 {
namespace elements {
namespace materialmodels {

/**
 * @brief Computes the stress by differentiation of the energy.
 */
template <class MaterialModel_> struct AutomaticStressTrait {
  template <class MaterialModel>
  typename MaterialModel::Stress
  operator()(const MaterialModel &model,
             const typename MaterialModel::size_type id,
             const typename MaterialModel::DisplacementGradient &gradient,
             const typename MaterialModel::Time time) const noexcept {
    using size_type = typename MaterialModel::size_type;

    typename MaterialModel::Stress result;

    auto modified_gradient = gradient;
    for (auto row = size_type{0}; row < MaterialModel::degrees_of_freedom();
         ++row) {
      for (auto col = size_type{0}; col < MaterialModel::dimension(); ++col) {
        result[row][col] = tensor::differentiate(
            [&](const typename MaterialModel::value_type t) {
              modified_gradient[row][col] = t;
              return compute_energy(model, id, modified_gradient, time);
            },
            gradient[row][col]);
        modified_gradient[row][col] = gradient[row][col];
      }
    }

    return result;
  }
};

} // namespace materialmodels
} // namespace elements
} // namespace ae108
