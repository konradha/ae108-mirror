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
 * @brief Sets the mesh of all target entities as periodic copies of the source
 * entity, using the affine transformation.
 *
 * The function supportes line (dim=1) and surface (dim=2) entities. Target and
 * surface entity need to have the same dimension.
 *
 * @param target_entities Target entities.
 * @param source_entity Source entity.
 * @param affine_transform Affine transformation matrix (4x4 matrix, by row)
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L1423
 */
void set_periodic(const std::vector<std::pair<int, int>> &target_entities,
                  const std::pair<int, int> &source_entity,
                  const std::array<double, 16> &affine_transform) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108