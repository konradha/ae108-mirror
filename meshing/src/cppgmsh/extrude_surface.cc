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

#include "ae108/meshing/cppgmsh/extrude_surface.h"
#include "ae108/meshing/cppgmsh/copy_entities.h"
#include "ae108/meshing/cppgmsh/extrude_entities.h"
#include "ae108/meshing/cppgmsh/get_centroid_of.h"
#include "ae108/meshing/cppgmsh/rotate_entities.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include "ae108/meshing/cppgmsh/translate_entities.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::pair<int, int>
extrude_surface(const int surface_tag, const std::array<double, 3> &from,
                const std::array<double, 3> &to,
                const std::array<double, 3> &surface_normal) noexcept {

  auto crosssection = copy_entities({{2, surface_tag}});

  const auto translation = [&] {
    std::array<double, 3> from_to;
    Eigen::Map<Eigen::Matrix<double, 3, 1>>(from_to.data()) =
        (Eigen::Map<const Eigen::Vector3d>(to.data()) -
         Eigen::Map<const Eigen::Vector3d>(from.data()));
    return from_to;
  }();

  std::array<double, 3> axis;
  Eigen::AngleAxisd angleAxis;
  angleAxis = Eigen::Quaterniond().setFromTwoVectors(
      Eigen::Map<const Eigen::Matrix<double, 3, 1>>(surface_normal.data()),
      Eigen::Map<const Eigen::Matrix<double, 3, 1>>(translation.data()));
  Eigen::Map<Eigen::Matrix<double, 3, 1>>(axis.data()) = angleAxis.axis();

  rotate_entities(crosssection, {0, 0, 0}, axis, angleAxis.angle());

  translate_entities(crosssection, from, false);

  const auto extruded_entities = extrude_entities(crosssection, translation);

  return extruded_entities[1];
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108