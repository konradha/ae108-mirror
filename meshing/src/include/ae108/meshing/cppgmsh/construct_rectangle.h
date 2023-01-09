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

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Constructs a rectangle, with lower left corner at {x, y, z} and
 * upper right corner at {x+dx,y+dy,z}. Returns the gmsh entity.
 *
 * @param origin Position of the rectangle origin {x,y,z}.
 * @param side_lengths Side lengths: {dx,dy}, {dy,dz} or {dx,dz}, according to
 * plane.
 * @param rounded_radius If > 0 the corners are rounded.
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L2229
 */
std::pair<int, int>
construct_rectangle(const std::array<double, 3> origin,
                    const std::array<double, 2> side_lengths,
                    const double rounded_radius = 0.) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108