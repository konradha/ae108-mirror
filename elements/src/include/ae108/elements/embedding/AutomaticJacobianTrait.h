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