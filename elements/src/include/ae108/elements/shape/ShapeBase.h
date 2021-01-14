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

#include "ae108/elements/tensor/Tensor.h"

namespace ae108 {
namespace elements {
namespace shape {

template <class SizeType_, class ValueType_, SizeType_ Dimension_,
          SizeType_ Size_>
class ShapeBase {
public:
  using size_type = SizeType_;
  using value_type = ValueType_;

  /**
   * @brief The number of shape functions.
   */
  static constexpr size_type size() noexcept { return Size_; }

  /**
   * @brief The dimension of reference space.
   */
  static constexpr size_type dimension() noexcept { return Dimension_; }

  /**
   * @brief Type of the gradient of a shape function.
   */
  using Gradient = tensor::Tensor<value_type, Dimension_>;

  /**
   * @brief Type of a point in reference space.
   */
  using Point = tensor::Tensor<value_type, Dimension_>;

  /**
   * @brief A collection of entities of size size().
   */
  template <class T> using Collection = tensor::Tensor<T, Size_>;

protected:
  ~ShapeBase() = default;
};

} // namespace shape
} // namespace elements
} // namespace ae108