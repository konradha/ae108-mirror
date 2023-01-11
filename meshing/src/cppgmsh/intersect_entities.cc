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

#include "ae108/meshing/cppgmsh/intersect_entities.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::vector<std::pair<int, int>>
intersect_entities(const std::vector<std::pair<int, int>> &object_entities,
                   const std::vector<std::pair<int, int>> &tool_entities,
                   const bool remove_object, const bool remove_tool) noexcept {
  std::vector<gmsh::vectorpair> temp;
  gmsh::vectorpair intersected_entities;
  gmsh::model::occ::intersect(object_entities, tool_entities,
                              intersected_entities, temp, -1, remove_object,
                              remove_tool);
  return intersected_entities;
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108