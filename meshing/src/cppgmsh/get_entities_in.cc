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

#include "ae108/meshing/cppgmsh/get_entities_in.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::vector<std::pair<int, int>>
get_entities_in(const BoundingBox<std::array<double, 3>> &bounding_box,
                const int dimension, const double tol) noexcept {

  gmsh::vectorpair entities;
  gmsh::model::getEntitiesInBoundingBox(
      bounding_box.min[0] - tol, bounding_box.min[1] - tol,
      bounding_box.min[2] - tol, bounding_box.max[0] + tol,
      bounding_box.max[1] + tol, bounding_box.max[2] + tol, entities,
      dimension);
  return entities;
}

std::vector<std::pair<int, int>>
get_entities_in(const std::pair<int, int> &physical_group) noexcept {

  std::vector<int> tags;
  gmsh::model::getEntitiesForPhysicalGroup(physical_group.first,
                                           physical_group.second, tags);

  std::vector<std::pair<int, int>> entities;
  entities.reserve(tags.size());
  for (const auto &tag : tags)
    entities.push_back({physical_group.first, tag});

  return entities;
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108