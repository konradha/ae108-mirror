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

#include "ae108/meshing/cppgmsh/get_elements_in.h"
#include <gmsh.h>
#include <vector>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::vector<std::pair<int, std::size_t>>
get_elements_in(const std::pair<int, int> &entity) noexcept {
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t>> elementTags;
  std::vector<std::vector<std::size_t>> nodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags,
                                 entity.first, entity.second);

  std::vector<std::pair<int, std::size_t>> elements;
  for (std::size_t i = 0; i < elementTypes.size(); i++)
    for (std::size_t j = 0; j < elementTags[i].size(); j++)
      elements.push_back({elementTypes[i], elementTags[i][j]});

  return elements;
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108