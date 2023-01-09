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

#include "ae108/meshing/cppgmsh/set_domain_entities_periodic.h"
#include "ae108/meshing/cppgmsh/as_affine_transform.h"
#include "ae108/meshing/cppgmsh/get_centroid_of.h"
#include "ae108/meshing/cppgmsh/get_entities_in.h"
#include "ae108/meshing/cppgmsh/set_periodic.h"
#include <Eigen/Dense>
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

void set_domain_entities_periodic(
    const BoundingBox<std::array<double, 3>> &source,
    const BoundingBox<std::array<double, 3>> &target, const int dim,
    const double tol) noexcept {

  assert(dim == 1 || dim == 2);

  std::array<double, 3> translation;
  Eigen::Map<Eigen::Vector3d>(translation.data()) =
      Eigen::Map<const Eigen::Vector3d>(target.min.data()) -
      Eigen::Map<const Eigen::Vector3d>(source.min.data());

  for (const auto &source_surface : get_entities_in(source, dim))
    for (const auto &target_surface : get_entities_in(target, dim))
      if (Eigen::Map<const Eigen::Vector3d>(
              get_centroid_of(target_surface).data())
              .isApprox(
                  Eigen::Map<const Eigen::Vector3d>(
                      get_centroid_of(source_surface).data()) +
                      Eigen::Map<const Eigen::Vector3d>(translation.data()),
                  tol))
        set_periodic({{2, target_surface.second}}, {2, source_surface.second},
                     as_affine_transform(translation));
}
} // namespace cppgmsh
} // namespace meshing
} // namespace ae108