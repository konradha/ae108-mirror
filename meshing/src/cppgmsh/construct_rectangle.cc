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

#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::pair<int, int>
construct_rectangle(const std::array<double, 3> origin,
                    const std::array<double, 2> side_lengths,
                    const double rounded_radius) noexcept {
  return {2, gmsh::model::occ::addRectangle(origin[0], origin[1], origin[2],
                                            side_lengths[0], side_lengths[1],
                                            -1, rounded_radius)};
}
} // namespace cppgmsh
} // namespace meshing
} // namespace ae108