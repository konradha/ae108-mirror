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
#include <cstddef>
#include <functional>

namespace ae108 {
namespace elements {
namespace integrator {

template <class SizeType_, class ValueType_, class RealType_, SizeType_ Size_,
          SizeType_ Dimension_>
struct IntegratorBase {
  using size_type = SizeType_;
  using value_type = ValueType_;
  using real_type = RealType_;

  /**
   * @brief Contains vectors of function evaluations at each selected point.
   */
  template <size_type Dim_>
  using DiscretizedFunction = tensor::Tensor<value_type, Size_, Dim_>;

  /**
   * @brief The integrand is evaluated at values of type Point, e.g. the
   * displacement gradient.
   */
  template <size_type Dim_>
  using Point = tensor::Tensor<value_type, Dim_, Dimension_>;

  /**
   * @brief The number of points a discretized function is evaluated at, e.g.
   * the number of shape functions.
   */
  static constexpr size_type size() noexcept { return Size_; }

  /**
   * @brief The integrand is evaluated at matrices Point. This method returns
   * the number of columns of these matrices. In the case of the displacement
   * gradient, this is the dimension of physical space.
   */
  static constexpr size_type dimension() noexcept { return Dimension_; }

protected:
  ~IntegratorBase() = default;
};

} // namespace integrator
} // namespace elements
} // namespace ae108