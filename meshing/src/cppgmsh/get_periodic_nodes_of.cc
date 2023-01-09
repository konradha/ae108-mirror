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

#include "ae108/meshing/cppgmsh/get_periodic_nodes_of.h"
#include <cassert>
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
get_periodic_nodes_of(const std::pair<int, int> &target_entity,
                      std::pair<int, int> *source_entity,
                      std::array<double, 16> *affine_transform) noexcept {

  assert(target_entity.first == 1 || target_entity.first == 2);

  std::vector<size_t> targetNodeTags, sourceNodeTags;
  std::vector<double> _;

  int source_tag;
  std::vector<double> transform;
  std::pair<std::vector<std::size_t>, std::vector<std::size_t>> node_tags;
  gmsh::model::mesh::getPeriodicNodes(target_entity.first, target_entity.second,
                                      source_tag, node_tags.first,
                                      node_tags.second, transform, true);
  if (source_entity) {
    source_entity->first = target_entity.first;
    source_entity->second = source_tag;
  }

  if (affine_transform)
    std::move(transform.data(), transform.data() + 16,
              affine_transform->begin());

  return node_tags;
};

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108