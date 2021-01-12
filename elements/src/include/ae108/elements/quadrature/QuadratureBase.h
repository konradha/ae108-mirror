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

#include <array>
#include <cstddef>

namespace ae108 {
namespace elements {
namespace quadrature {

template <class SizeType_, class PointType_, class WeightType_, SizeType_ Size_>
struct QuadratureBase {
  using size_type = SizeType_;

  /**
   * @brief The type of the quadrature points.
   */
  using Point = PointType_;

  /**
   * @brief The type of the quadrature weights.
   */
  using Weight = WeightType_;

  /**
   * @brief A collection of entities of size size().
   */
  template <class T> using Collection = std::array<T, Size_>;

  /**
   * @brief The number of quadrature points.
   */
  static constexpr size_type size() noexcept { return Size_; }

  struct Data {
    Collection<Point> points;
    Collection<Weight> weights;
  };

protected:
  ~QuadratureBase() = default;
};

} // namespace quadrature
} // namespace elements
} // namespace ae108