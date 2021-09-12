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