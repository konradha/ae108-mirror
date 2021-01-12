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