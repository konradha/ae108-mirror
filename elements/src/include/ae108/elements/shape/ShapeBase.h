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