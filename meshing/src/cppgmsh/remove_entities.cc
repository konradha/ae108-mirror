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

#include "ae108/meshing/cppgmsh/remove_entities.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

void remove_entities(const std::vector<std::pair<int, int>> &entities,
                     const bool recursive) noexcept {
  gmsh::model::occ::remove(entities, recursive);
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108