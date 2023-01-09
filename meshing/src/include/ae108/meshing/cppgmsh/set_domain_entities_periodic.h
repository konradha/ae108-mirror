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
 * @brief Set the meshes of the entities in domain target as
 * periodic copies of the meshes of entities in domain source
 *
 * @param source Source domain.
 * @param target Target domain.
 * @param dim Dimension of the entities (1 or 2).
 * @param tol Spatial tolerance for entity matching.
 */
void set_domain_entities_periodic(
    const BoundingBox<std::array<double, 3>> &source,
    const BoundingBox<std::array<double, 3>> &target, const int dim,
    const double tol = 1e-9) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108