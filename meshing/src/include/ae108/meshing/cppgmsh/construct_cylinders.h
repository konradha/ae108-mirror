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
#include <vector>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Constructs cylinders of differing radius from each segment using gmsh.
 * Returns the entity.
 *
 * @param nodes Node positions.
 * @param connectivity Nodal connectivity of each segment.
 * @param radii Radii of all segments.
 * @param capped Option to add a cap on each cylinder end.
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L2414
 */
std::vector<std::pair<int, int>>
construct_cylinders(const std::vector<std::array<double, 3>> &nodes,
                    const std::vector<std::array<std::size_t, 2>> &connectivity,
                    const std::vector<double> &radii,
                    const bool capped = true) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108