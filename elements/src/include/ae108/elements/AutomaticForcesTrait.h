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

#include "ae108/elements/compute_energy.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include "ae108/elements/tensor/differentiate.h"

namespace ae108 {
namespace elements {

/**
 * @brief Computes the forces by differentiation of the energy.
 */
template <class Element_> struct AutomaticForcesTrait {
  template <class Element>
  typename Element::Forces
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &time) const noexcept {
    auto result = typename Element::Forces();
    auto result_matrix = tensor::as_matrix_of_rows(&result);

    auto modified_u = u;
    for (auto row = Eigen::Index{0}; row < result_matrix.rows(); ++row) {
      for (auto col = Eigen::Index{0}; col < result_matrix.cols(); ++col) {
        result_matrix(row, col) = tensor::differentiate(
            [&](const typename Element::value_type t) {
              modified_u[row][col] = t;
              return compute_energy(element, modified_u, time);
            },
            u[row][col]);
        modified_u[row][col] = u[row][col];
      }
    }

    return result;
  }
};

} // namespace elements
} // namespace ae108