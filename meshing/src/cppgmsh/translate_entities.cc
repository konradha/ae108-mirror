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

#include "ae108/meshing/cppgmsh/translate_entities.h"
#include "ae108/meshing/cppgmsh/copy_entities.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::vector<std::pair<int, int>>
translate_entities(std::vector<std::pair<int, int>> entities,
                   const std::array<double, 3> &translation,
                   const bool copy) noexcept {
  if (copy)
    entities = copy_entities(entities);

  gmsh::model::occ::translate(entities, translation[0], translation[1],
                              translation[2]);

  return entities;
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108