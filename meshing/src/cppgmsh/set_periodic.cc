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

#include "ae108/meshing/cppgmsh/set_periodic.h"
#include <cassert>
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

void set_periodic(const std::vector<std::pair<int, int>> &target_entities,
                  const std::pair<int, int> &source_entity,
                  const std::array<double, 16> &affine_transform) noexcept {

  assert(source_entity.first == 1 || source_entity.first == 2);

  std::vector<int> target_tags;
  for (const auto &target_entity : target_entities) {
    assert(source_entity.first == target_entity.first);
    target_tags.push_back(target_entity.second);
  }

  std::vector<double> transform(affine_transform.begin(),
                                affine_transform.end());

  gmsh::model::mesh::setPeriodic(source_entity.first, target_tags,
                                 {source_entity.second}, transform);
};

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108