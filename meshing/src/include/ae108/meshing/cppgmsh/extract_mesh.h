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
#include <map>
#include <vector>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Extract mesh from Gmsh.
 *
 * @param element_dimension dimension of the elements of interest.
 */
template <std::size_t coordinate_dimension>

std::tuple<std::vector<std::array<double, coordinate_dimension>>,
           std::vector<std::vector<std::size_t>>,
           std::map<std::size_t, std::size_t>,
           std::map<std::size_t, std::size_t>>
extract_mesh(
    const std::size_t element_dimension = coordinate_dimension) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108