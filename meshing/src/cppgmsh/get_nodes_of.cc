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

#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

template <std::size_t Dimension>
std::vector<Node<Dimension>>
get_nodes_of(const std::pair<int, int> &entity,
             const bool remove_duplicates) noexcept {

  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> _;
  gmsh::model::mesh::getNodes(nodeTags, coord, _, entity.first, entity.second,
                              true, false);

  std::vector<Node<Dimension>> nodes(nodeTags.size());
  for (std::size_t i = 0; i < nodeTags.size(); i++) {
    nodes[i].id = nodeTags[i];
    std::copy_n(coord.begin() + i * 3, Dimension, nodes[i].position.begin());
  }

  if (remove_duplicates) {
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  }

  return nodes;
}

template std::vector<Node<3>>
get_nodes_of<3>(const std::pair<int, int> &entity,
                const bool keep_duplicates) noexcept;

template std::vector<Node<2>>
get_nodes_of<2>(const std::pair<int, int> &entity,
                const bool keep_duplicates) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108