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

#include "ae108/meshing/cppgmsh/get_centroid_of.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::array<double, 3>
get_centroid_of(const std::pair<int, int> &entity) noexcept {
  std::array<double, 3> centroid;
  gmsh::model::occ::getCenterOfMass(entity.first, entity.second, centroid[0],
                                    centroid[1], centroid[2]);
  return centroid;
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108