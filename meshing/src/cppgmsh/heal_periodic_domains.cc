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

#include "ae108/meshing/cppgmsh/heal_periodic_domains.h"
#include "ae108/meshing/cppgmsh/get_centroid_of.h"
#include "ae108/meshing/cppgmsh/get_entities_in.h"
#include "ae108/meshing/cppgmsh/heal_periodic_surfaces.h"
#include <Eigen/Dense>
#include <array>
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

void heal_periodic_domains(const BoundingBox<std::array<double, 3>> &source,
                           const BoundingBox<std::array<double, 3>> &target,
                           const double tol) noexcept {
  std::array<double, 3> translation;
  Eigen::Map<Eigen::Vector3d>(translation.data()) =
      Eigen::Map<const Eigen::Vector3d>(target.min.data()) -
      Eigen::Map<const Eigen::Vector3d>(source.min.data());
  for (const auto &source_surface : get_entities_in(source, 2)) {
    const auto source_centroid = get_centroid_of(source_surface);
    for (const auto &target_surface : get_entities_in(target, 2)) {
      const auto target_centroid = get_centroid_of(target_surface);
      if (Eigen::Map<const Eigen::Vector3d>(target_centroid.data())
              .isApprox(
                  Eigen::Map<const Eigen::Vector3d>(source_centroid.data()) +
                      Eigen::Map<const Eigen::Vector3d>(translation.data()),
                  tol))
        heal_periodic_surfaces(source_surface.second, target_surface.second,
                               translation);
    }
  }
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108