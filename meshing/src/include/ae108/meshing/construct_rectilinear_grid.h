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

/**
 * @brief Constructs a d-dimensional rectilinear grid based on the (non-)uniform
 * input coordinates. Returns points and edges.
 *
 * @param coordinates d sets of coordinates that span the grid.
 * @note https://en.wikipedia.org/wiki/Regular_grid
 */

template <std::size_t Dimension>
std::tuple<std::vector<std::array<double, Dimension>>,
           std::vector<std::array<std::size_t, 2>>>
construct_rectilinear_grid(
    const std::array<std::vector<double>, Dimension> &coordinates) noexcept;

} // namespace meshing
} // namespace ae108