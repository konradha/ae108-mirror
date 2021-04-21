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