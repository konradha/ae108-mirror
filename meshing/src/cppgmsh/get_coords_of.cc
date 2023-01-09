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

#include "ae108/meshing/cppgmsh/get_coords_of.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::array<double, 3> get_coords_of(const int point_tag) noexcept {
  std::vector<double> _, coord;
  gmsh::model::getValue(0, point_tag, _, coord);
  return {coord[0], coord[1], coord[2]};
}

std::vector<std::array<double, 3>>
get_coords_of(const std::vector<int> &point_tags) noexcept {
  std::vector<std::array<double, 3>> coords(point_tags.size());
  for (std::size_t i = 0; i < point_tags.size(); i++)
    coords[i] = get_coords_of(point_tags[i]);
  return coords;
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108