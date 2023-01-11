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

#include "ae108/meshing/cppgmsh/Node.h"
#include <vector>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Get the nodes of one or more entities. Optionally avoids
 * duplicate nodes (e.g. on entity boundaries). Returns a vector of Nodes.
 *
 * @param dim The dimension of the model entities of interest. Set to < 0 to get
 * all nodes of the mesh.
 * @param tag The tag of the model entity of interest. Set to < 0 to get the
 * nodes of all entities of dimension dim.
 * @param remove_duplicates If true, node duplicates are removed.
 * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L662
 */

template <std::size_t Dimension>
std::vector<Node<Dimension>>
get_nodes_of(const std::pair<int, int> &entity,
             const bool remove_duplicates = true) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108