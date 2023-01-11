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

#include <vector>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Returns vector of intersected entities.
 *
 * @param object_entities Object entities.
 * @param tool_entities Tool entities.
 * @param remove_object If true, the object entities will be removed after the
 * intersection was performed.
 * @param remove_tool If true, the tool entities will be removed after the
 * intersection was performed.
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L2616
 */
std::vector<std::pair<int, int>>
intersect_entities(const std::vector<std::pair<int, int>> &object_entities,
                   const std::vector<std::pair<int, int>> &tool_entities,
                   const bool remove_object = true,
                   const bool remove_tool = true) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108