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