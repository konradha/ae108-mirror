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

#include "ae108/elements/materialmodels/compute_stress.h"
#include "ae108/elements/tensor/differentiate.h"

namespace ae108 {
namespace elements {
namespace materialmodels {

/**
 * @brief Computes the forces by differentiation of the stresses.
 */
template <class MaterialModel_> struct AutomaticTangentMatrixTrait {
  template <class MaterialModel>
  typename MaterialModel::TangentMatrix
  operator()(const MaterialModel &model,
             const typename MaterialModel::size_type id,
             const typename MaterialModel::DisplacementGradient &gradient,
             const typename MaterialModel::Time time) const noexcept {
    using size_type = typename MaterialModel::size_type;

    typename MaterialModel::TangentMatrix result;

    auto modified_gradient = gradient;
    for (auto i = size_type{0}; i < MaterialModel::degrees_of_freedom(); ++i) {
      for (auto j = size_type{0}; j < MaterialModel::dimension(); ++j) {
        for (auto k = size_type{0}; k < MaterialModel::degrees_of_freedom();
             ++k) {
          for (auto l = size_type{0}; l < MaterialModel::dimension(); ++l) {
            result[i][j][k][l] = tensor::differentiate(
                [&](const typename MaterialModel::value_type t) {
                  modified_gradient[k][l] = t;
                  return compute_stress(model, id, modified_gradient,
                                        time)[i][j];
                },
                gradient[k][l]);
            modified_gradient[k][l] = gradient[k][l];
          }
        }
      }
    }

    return result;
  }
};

} // namespace materialmodels
} // namespace elements
} // namespace ae108
