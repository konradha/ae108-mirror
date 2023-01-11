// © 2022 ETH Zurich, Mechanics and Materials Lab
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
 * @brief Returns boundary entities (dim-1) of the specified entities, or
 * optionally the boundary points, if recursive = true.
 *
 * @param entities Entities of interest.
 * @param recursive If true, returns boundary points (dim=0).
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L284
 */
std::vector<std::pair<int, int>>
get_boundary_of(const std::vector<std::pair<int, int>> &entities,
                const bool recursive = false) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108