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

namespace ae108 {
namespace elements {
namespace shape {

/**
 * @brief Contains a member value that defines the points of the shape.
 */
template <class Shape> struct PointTrait;

#define AE108_ELEMENTS_SHAPE_DEFINE_POINTS(name, ...)                          \
  template <> struct PointTrait<name> {                                        \
    const typename name::template Collection<typename name::Point> &           \
    operator()(const name &) noexcept {                                        \
      static constexpr                                                         \
          typename name::template Collection<typename name::Point>             \
              points = __VA_ARGS__;                                            \
      return points;                                                           \
    }                                                                          \
  }

} // namespace shape
} // namespace elements
} // namespace ae108