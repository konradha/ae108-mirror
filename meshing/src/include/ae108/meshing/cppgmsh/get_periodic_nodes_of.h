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
 * @brief Returns the corresponding node tags of the target entity (out.first)
 * and the source entity (out.second)
 *
 * @param target_entity Target entity.
 * @param source_entity Optional output of the source entity.
 * @param affine_transform Optional output of the affine transformation.
 */
std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
get_periodic_nodes_of(
    const std::pair<int, int> &target_entity,
    std::pair<int, int> *source_entity = nullptr,
    std::array<double, 16> *affine_transform = nullptr) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108