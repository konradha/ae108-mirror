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

#include "ae108/meshing/construct_voronoi_cell.h"
#include <Eigen/Dense>
#include <range/v3/action/sort.hpp>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/enumerate.hpp>
#include <voro++.hh>

namespace ae108 {
namespace meshing {

template <>
BoundaryRepresentation<std::size_t, double, 3>
construct_voronoi_cell<std::size_t, double, 3>(
    const std::vector<std::array<double, 3>> &point_cloud,
    std::vector<std::pair<std::size_t, std::size_t>> *periodic_faces) noexcept {

  assert(point_cloud.size());

  voro::voronoicell voronoi_cell;
  const auto bound = 1e3;
  voronoi_cell.init(-bound, bound, -bound, bound, -bound, bound);
  for (auto &point : point_cloud)
    voronoi_cell.plane(point[0], point[1], point[2]);

  const auto vertices = [](voro::voronoicell &voronoi_cell) {
    std::vector<double> vertices;
    voronoi_cell.vertices(vertices);
    return vertices;
  }(voronoi_cell);

  const auto faces = [](voro::voronoicell &voronoi_cell) {
    std::vector<std::vector<std::array<std::size_t, 2>>> faces;
    std::vector<int> temp;
    voronoi_cell.face_vertices(temp);
    for (int i = 0; i < (int)temp.size(); i = i + temp[i] + 1) {
      std::vector<std::array<std::size_t, 2>> edges;
      for (int j = 1; j < temp[i]; j++) {
        edges.push_back(
            {std::size_t(temp[i + j]), std::size_t(temp[i + j + 1])});
        ranges::sort(edges.back());
      }
      edges.push_back(
          {std::size_t(temp[i + temp[i]]), std::size_t(temp[i + 1])});
      ranges::sort(edges.back());
      faces.push_back(edges);
    }
    return faces;
  }(voronoi_cell);

  const auto face_normals = [](voro::voronoicell &voronoi_cell) {
    std::vector<BoundaryRepresentation<std::size_t, double, 3>::Point>
        face_normals(voronoi_cell.number_of_faces());
    std::vector<double> temp;
    voronoi_cell.normals(temp);
    for (const auto &[index, normal] :
         ranges::views::enumerate(temp | ranges::views::chunk(3)))
      ranges::copy(normal, face_normals[index].begin());
    return face_normals;
  }(voronoi_cell);

  auto add_edge =
      [](const BoundaryRepresentation<std::size_t, double, 3>::Edge &edge,
         BoundaryRepresentation<std::size_t, double, 3>::Edges &edges)
      -> std::size_t {
    for (auto [index, connectivity] : ranges::views::enumerate(edges))
      if (edge == connectivity)
        return index;
    edges.push_back(edge);
    return edges.size() - 1;
  };

  BoundaryRepresentation<std::size_t, double, 3> brep;
  brep.faces.resize(voronoi_cell.number_of_faces());
  for (const auto &[index, face] : ranges::views::enumerate(faces))
    for (auto &edge : face)
      brep.faces[index].push_back(add_edge(edge, brep.edges));

  brep.vertices.resize(vertices.size() / 3);
  for (const auto &[index, vertex] :
       ranges::views::enumerate(vertices | ranges::views::chunk(3)))
    ranges::copy(vertex, brep.vertices[index].begin());

  if (periodic_faces)
    for (std::size_t i = 0; i < face_normals.size(); i++)
      for (std::size_t j = i + 1; j < face_normals.size(); j++)
        if (Eigen::Map<const Eigen::Vector3d>(face_normals[i].data())
                .cross(
                    Eigen::Map<const Eigen::Vector3d>(face_normals[j].data()))
                .isZero())
          periodic_faces->push_back({i, j});

  return brep;
}

} // namespace meshing
} // namespace ae108