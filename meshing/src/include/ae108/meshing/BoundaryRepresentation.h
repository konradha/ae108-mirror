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
#include <array>
#include <set>
#include <vector>

namespace ae108 {
namespace meshing {
template <class SizeType_, class ValueType_, SizeType_ Dimension_>
struct BoundaryRepresentation {
  using Point = std::array<ValueType_, Dimension_>;

  /**
   * @brief The position per vertex.
   */
  using Vertices = std::vector<Point>;

  /**
   * @brief Indicies of two vertices that form an edge.
   */
  using Edge = std::array<SizeType_, 2>;

  /**
   * @brief Vector of all edges.
   */
  using Edges = std::vector<Edge>;

  /**
   * @brief Edges that enclose a face.
   */
  using Faces = std::vector<std::vector<SizeType_>>;

  /**
   * @brief Returns the dimension.
   */
  static constexpr SizeType_ dimension() noexcept { return Dimension_; }

  /**
   * @brief Returns all vertices of an edge.
   */
  std::vector<Point> vertices_of_edge(SizeType_ edge) const noexcept {
    std::vector<Point> edge_vertices;
    for (const auto &vertex : edges[edge])
      edge_vertices.push_back(vertices[vertex]);
    return edge_vertices;
  };

  /**
   * @brief Returns all vertices of a face.
   */
  std::vector<Point> vertices_of_face(SizeType_ index) const noexcept {
    std::set<SizeType_> face_vertex_indices;
    for (const auto &edge : faces[index])
      for (const auto &vertex : edges[edge])
        face_vertex_indices.insert(vertex);

    std::vector<Point> face_vertices;
    face_vertices.reserve(face_vertex_indices.size());
    for (const auto &vertex_index : face_vertex_indices)
      face_vertices.push_back(vertices[vertex_index]);

    return face_vertices;
  };

  Vertices vertices;
  Edges edges;
  Faces faces;
};
} // namespace meshing
} // namespace ae108