// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/meshing/BoundingBox.h"
#include <cassert>
#include <vector>

namespace ae108 {
namespace meshing {

/**
 * @brief Returns the bounding box of points
 *
 * @param points A vector of points.
 */

template <class Point>
BoundingBox<Point> bounding_box_of(const std::vector<Point> &points) noexcept {
  assert(points.size());
  Point min(points[0]);
  Point max(points[0]);
  for (const auto &point : points)
    for (std::size_t i = 0; i < point.size(); i++) {
      min[i] = std::min(point[i], min[i]);
      max[i] = std::max(point[i], max[i]);
    }
  return {min, max};
}

} // namespace meshing
} // namespace ae108