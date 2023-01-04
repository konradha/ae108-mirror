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
