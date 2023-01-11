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

#include "ae108/meshing/cppgmsh/heal_periodic_surfaces.h"
#include "ae108/meshing/cppgmsh/get_coords_of.h"
#include "ae108/meshing/cppgmsh/get_points_of.h"
#include <Eigen/Dense>
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

void heal_periodic_surfaces(const int source_surface, const int target_surface,
                            const std::array<double, 3> &translation,
                            const double tol) noexcept {

  for (const auto &source : get_points_of({2, source_surface})) {

    std::array<double, 3> query;
    Eigen::Map<Eigen::Vector3d>(query.data()) =
        (Eigen::Map<const Eigen::Vector3d>(get_coords_of(source).data()) +
         Eigen::Map<const Eigen::Vector3d>(translation.data()))
            .eval();

    bool target_found = false;
    for (const auto &target : get_points_of({2, target_surface}))
      target_found +=
          Eigen::Map<const Eigen::Vector3d>(get_coords_of(target).data())
              .isApprox(Eigen::Map<const Eigen::Vector3d>(query.data()), tol);

    if (!target_found) {
      const auto missing_point =
          gmsh::model::occ::addPoint(query[0], query[1], query[2]);

      std::vector<int> volumes, lines;
      gmsh::model::getAdjacencies(2, target_surface, volumes, lines);

      gmsh::vectorpair ov;
      std::vector<gmsh::vectorpair> ovv;
      if (volumes.size())
        gmsh::model::occ::fragment(
            [](const std::vector<int> &volumes) {
              std::vector<std::pair<int, int>> volume_entities;
              for (const auto &tag : volumes)
                volume_entities.push_back({3, tag});
              return volume_entities;
            }(volumes),
            {{0, missing_point}}, ov, ovv);
      else
        gmsh::model::occ::fragment({{2, target_surface}}, {{0, missing_point}},
                                   ov, ovv);
    }
  }
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108