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

#include "ae108/meshing/cppgmsh/construct_polyhedron.h"
#include <gmsh.h>
#include <map>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

template <>
std::pair<int, int> construct_polyhedron(
    const BoundaryRepresentation<std::size_t, double, 3> &brep) noexcept {

  std::map<std::size_t, int> point2tag;
  for (auto &&[index, vertex] : ranges::views::enumerate(brep.vertices))
    point2tag[index] =
        gmsh::model::occ::addPoint(vertex[0], vertex[1], vertex[2]);

  std::map<std::size_t, int> edge2tag;
  for (auto &&[index, edge] : ranges::views::enumerate(brep.edges))
    edge2tag[index] =
        gmsh::model::occ::addLine(point2tag.at(edge[0]), point2tag.at(edge[1]));

  std::vector<int> surfaceTags;
  for (const auto &line_loop : brep.faces)
    surfaceTags.push_back(
        gmsh::model::occ::addPlaneSurface({gmsh::model::occ::addCurveLoop(
            line_loop | ranges::view::transform([&edge2tag](int edge) {
              return edge2tag.at(edge);
            }) |
            ranges::to<std::vector>)}));

  return {3, gmsh::model::occ::addVolume(
                 {gmsh::model::occ::addSurfaceLoop(surfaceTags)})};
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108