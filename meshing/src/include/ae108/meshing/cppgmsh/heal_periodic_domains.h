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

#include "ae108/meshing/BoundingBox.h"
#include <array>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Matches points between all periodic surfaces in the source and target
 * domain. If a point on the source surface is missing on the target surface,
 * it will be added to the target surface. Assumption: Source and target domain
 * are translational periodic.
 *
 * @param source Source domain spanned by a bounding box.
 * @param target Target domain spanned by a bounding box.
 * @param tol Tolerance for point matching.
 */
void heal_periodic_domains(const BoundingBox<std::array<double, 3>> &source,
                           const BoundingBox<std::array<double, 3>> &target,
                           const double tol = 1e-9) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108