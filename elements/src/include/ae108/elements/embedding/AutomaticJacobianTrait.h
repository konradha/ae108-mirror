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

#include "ae108/elements/embedding/embed_point.h"
#include "ae108/elements/tensor/Tensor.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include "ae108/elements/tensor/differentiate.h"

namespace ae108 {
namespace elements {
namespace embedding {

template <class Embedding_> struct AutomaticJacobianTrait {
  template <class Embedding>
  using Jacobian = tensor::Tensor<typename Embedding::value_type,
                                  Embedding::physical_dimension(),
                                  Embedding::reference_dimension()>;

  template <class Embedding>
  Jacobian<Embedding>
  operator()(const Embedding &embedding,
             const typename Embedding::ReferencePoint &point) const noexcept {
    auto result = Jacobian<Embedding>();
    auto result_matrix = tensor::as_matrix_of_rows(&result);

    for (auto row = Eigen::Index{0}; row < result_matrix.rows(); ++row) {
      for (auto col = Eigen::Index{0}; col < result_matrix.cols(); ++col) {
        auto modified_point = point;
        result_matrix(row, col) = tensor::differentiate(
            [&](const typename Embedding::value_type t) {
              modified_point[col] = t;
              return embed_point<Embedding>(embedding, modified_point)[row];
            },
            point[col]);
      }
    }

    return result;
  }
};

} // namespace embedding
} // namespace elements
} // namespace ae108