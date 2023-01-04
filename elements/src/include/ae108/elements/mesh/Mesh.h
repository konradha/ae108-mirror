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

#include <vector>

namespace ae108 {
namespace elements {
namespace mesh {

template <class Point_> class Mesh {
public:
  using Point = Point_;
  using size_type = typename std::vector<Point>::size_type;

  /**
   * @brief The vertex indices per element.
   */
  using Connectivity = std::vector<std::vector<size_type>>;

  /**
   * @brief The position per vertex.
   */
  using Positions = std::vector<Point>;

  /**
   * @brief Create a mesh from a connectivity and positions.
   * @param positions Should contain a position for every vertex referenced in
   * the connectivity.
   */
  explicit Mesh(Connectivity connectivity, Positions positions) noexcept
      : connectivity_(std::move(connectivity)),
        positions_(std::move(positions)) {}

  /**
   * @brief Returns the connectivity of the mesh.
   */
  const Connectivity &connectivity() const noexcept { return connectivity_; }

  /**
   * @brief Returns the position of the vertex with the given index.
   * @throw std::out_of_range If the mesh does not contain position information
   * for the given index.
   */
  const Point &position_of_vertex(const size_type index) const {
    return positions_.at(index);
  }

  /**
   * @brief Returns the number of positions in this mesh.
   */
  size_type number_of_positions() const noexcept { return positions_.size(); }

private:
  Connectivity connectivity_;
  Positions positions_;
};

} // namespace mesh
} // namespace elements
} // namespace ae108