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
 * @brief Matches points between source and target surface. If a point on the
 * source surface is missing on the target surface, it will be added to the
 * target surface. Source and target surface are translational periodic.
 *
 * @param source Source surface tag.
 * @param target Target surface tag.
 * @param translation Translation vector from source to target.
 * @param tol Tolerance for point matching.
 *
 */
void heal_periodic_surfaces(const int source, const int target,
                            const std::array<double, 3> &translation,
                            const double tol = 1e-9) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108