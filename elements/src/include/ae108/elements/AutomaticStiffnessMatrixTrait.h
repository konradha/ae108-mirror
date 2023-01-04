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

#include "ae108/elements/compute_forces.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include "ae108/elements/tensor/differentiate.h"

namespace ae108 {
namespace elements {

/**
 * @brief Computes the stiffness matrix by differentiation of the forces.
 */
template <class Element_> struct AutomaticStiffnessMatrixTrait {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &time) const noexcept {
    using size_type = typename Element::size_type;
    constexpr auto size = Element::size();
    constexpr auto dofs = Element::degrees_of_freedom();

    auto result = typename Element::StiffnessMatrix();

    auto modified_u = u;
    for (auto node_row = size_type{0}; node_row < size; ++node_row) {
      for (auto dof_row = size_type{0}; dof_row < dofs; ++dof_row) {
        for (auto node_col = size_type{0}; node_col < size; ++node_col) {
          for (auto dof_col = size_type{0}; dof_col < dofs; ++dof_col) {
            result(node_row * dofs + dof_row, node_col * dofs + dof_col) =
                tensor::differentiate(
                    [&](const typename Element::value_type t) {
                      modified_u[node_col][dof_col] = t;
                      return compute_forces(element, modified_u,
                                            time)[node_row][dof_row];
                    },
                    u[node_col][dof_col]);
            modified_u[node_col][dof_col] = u[node_col][dof_col];
          }
        }
      }
    }

    return result;
  }
};

} // namespace elements
} // namespace ae108