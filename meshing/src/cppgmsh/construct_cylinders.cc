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

#include "ae108/meshing/cppgmsh/construct_cylinders.h"
#include <gmsh.h>
#include <range/v3/view/zip.hpp>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::pair<int, int> construct_cylinder(const std::array<double, 3> &p0,
                                       const std::array<double, 3> &p1,
                                       const double radius) {
  return {3,
          gmsh::model::occ::addCylinder(p0[0], p0[1], p0[2], p1[0] - p0[0],
                                        p1[1] - p0[1], p1[2] - p0[2], radius)};
}

std::pair<int, int> capped_cylinder(const std::array<double, 3> &p0,
                                    const std::array<double, 3> &p1,
                                    const double radius) {
  std::vector<gmsh::vectorpair> temp;
  gmsh::vectorpair cappedCylinder;
  gmsh::model::occ::fuse(
      {construct_cylinder(p0, p1, radius)},
      {{3, gmsh::model::occ::addSphere(p0[0], p0[1], p0[2], radius)},
       {3, gmsh::model::occ::addSphere(p1[0], p1[1], p1[2], radius)}},
      cappedCylinder, temp);
  return cappedCylinder[0];
}

std::vector<std::pair<int, int>>
construct_cylinders(const std::vector<std::array<double, 3>> &nodes,
                    const std::vector<std::array<std::size_t, 2>> &connectivity,
                    const std::vector<double> &radii,
                    const bool capped) noexcept {

  assert(connectivity.size() == radii.size());

  gmsh::vectorpair cylinder_entities;
  cylinder_entities.reserve(connectivity.size());
  for (const auto &&[segment, radius] :
       ranges::views::zip(connectivity, radii)) {

    const auto &p0 = nodes[segment[0]];
    const auto &p1 = nodes[segment[1]];

    if (capped)
      cylinder_entities.push_back(capped_cylinder(p0, p1, radius));
    else
      cylinder_entities.push_back(construct_cylinder(p0, p1, radius));
  }

  return cylinder_entities;
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108