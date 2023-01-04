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