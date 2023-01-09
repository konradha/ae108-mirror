// © 2022 ETH Zurich, Mechanics and Materials Lab
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
 * @brief Extrudes the given entities along a translation vector. Returns the
 * extruded entities.
 * @param entities Entities to extrude_entities.
 * @param translation Extrusion vector.
 * @param extrude_mesh If true, the mesh of the base area of the extruded bodies
 * will be extruded as well, i.e. the same layer of 3D elements will be repeated
 * several times. Otherwise, The bodies are meshed as general 3D bodies.
 * @param number_of_layers Number of element layers (if extrude_mesh is true).
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L2514
 */
std::vector<std::pair<int, int>>
extrude_entities(const std::vector<std::pair<int, int>> &entities,
                 const std::array<double, 3> &translation,
                 const bool extrude_mesh = false,
                 const int number_of_layers = 1) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108