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
 * @brief Returns affine transform (16 entries of a 4x4 matrix, by row) based on
 * translation.
 *
 * @param translation The 3D translation vector.
 */
std::array<double, 16>
as_affine_transform(const std::array<double, 3> &translation) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108