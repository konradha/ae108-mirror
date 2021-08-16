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
