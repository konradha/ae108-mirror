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
 * @brief Extrudes the surface entity from origin to target.
 *
 * @param surface_tag Tag of the surface defining the cross-section in xy-plane.
 * @param from Position of the origin.
 * @param to Position of the target.
 * @param surface_normal Normal of the surface.
 */
std::pair<int, int> extrude_surface(
    const int surface_tag, const std::array<double, 3> &from,
    const std::array<double, 3> &to,
    const std::array<double, 3> &surface_normal = {{0, 0, 1}}) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108