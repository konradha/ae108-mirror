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
#include <vector>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Returns all entities of dimension inside of a bounding box
 *
 * @param bounding_box The bounding box or search domain.
 * @param dimension Entity dimension, i.e 3-volume, 2-surface, 1-line, 0-point
 * @param tol Spatial tolerance.
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L307
 */
std::vector<std::pair<int, int>>
get_entities_in(const BoundingBox<std::array<double, 3>> &bounding_box,
                const int dimension, const double tol = 1e-6) noexcept;

/**
 * @brief Returns all entities in the specified physical group.
 *
 * @param physical_group The physical group {dim, tag}.
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L230
 */
std::vector<std::pair<int, int>>
get_entities_in(const std::pair<int, int> &physical_group) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108