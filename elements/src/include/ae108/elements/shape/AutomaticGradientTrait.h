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

#include "ae108/elements/shape/compute_values.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include "ae108/elements/tensor/differentiate.h"

namespace ae108 {
namespace elements {
namespace shape {

template <class Shape_> struct AutomaticGradientTrait {
  template <class Shape>
  typename Shape::template Collection<typename Shape::Gradient>
  operator()(const Shape &, const typename Shape::Point &xi) const noexcept {
    auto result =
        typename Shape::template Collection<typename Shape::Gradient>();
    auto result_matrix = tensor::as_matrix_of_rows(&result);

    for (auto row = Eigen::Index{0}; row < result_matrix.rows(); ++row) {
      for (auto col = Eigen::Index{0}; col < result_matrix.cols(); ++col) {
        auto modified_xi = xi;
        result_matrix(row, col) = tensor::differentiate(
            [&](const typename Shape::value_type t) {
              modified_xi[col] = t;
              return compute_values<Shape>(modified_xi)[row];
            },
            xi[col]);
      }
    }

    return result;
  }
};

} // namespace shape
} // namespace elements
} // namespace ae108