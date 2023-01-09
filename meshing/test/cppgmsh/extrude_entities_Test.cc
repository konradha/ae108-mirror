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

#include "ae108/meshing/cppgmsh/Context.h"
#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include "ae108/meshing/cppgmsh/extrude_entities.h"
#include "ae108/meshing/cppgmsh/fragment_entities.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_coords_of.h"
#include "ae108/meshing/cppgmsh/get_entities_of.h"
#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"

#include <cmath>
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::DoubleEq;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct extrude_Test : Test {};

TEST_F(extrude_Test, extrudes_rectangle) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle = construct_rectangle({0., 0., 0.}, {1., 1.});
  const std::array<double, 3> dr = {0.1, 0.2, 0.3};
  extrude_entities({rectangle}, dr);
  synchronize();

  const auto points = get_entities_of(0);
  ASSERT_THAT(points.size(), Eq(8));
  const auto edges = get_entities_of(1);
  ASSERT_THAT(edges.size(), Eq(12));
  const auto surfaces = get_entities_of(2);
  ASSERT_THAT(surfaces.size(), Eq(6));
  const auto volumes = get_entities_of(3);
  ASSERT_THAT(volumes.size(), Eq(1));

  for (std::size_t i = 0; i < 4; i++) {
    const auto coords_0 = get_coords_of(points[i].second);
    const auto coords_1 = get_coords_of(points[i + 4].second);
    for (std::size_t d = 0; d < 3; d++) {
      ASSERT_THAT(coords_0[d] + dr[d], DoubleEq(coords_1[d]));
    }
  }
}

TEST_F(extrude_Test, extrudes_mesh) {
  const auto gmshContext = Context(0, 0);
  const auto rectangle_1 = construct_rectangle({0., 0., 0.}, {1., 1.});
  const auto rectangle_2 = construct_rectangle({2., 0., 0.}, {1., 1.});
  const std::array<double, 3> dr = {0, 0, 1};
  auto extruded_entities_1 = extrude_entities({rectangle_1}, dr);
  std::pair<int, int> box_1 = {-1, -1};
  for (std::size_t i = 0; i < extruded_entities_1.size(); i++) {
    if (extruded_entities_1[i].first == 3) {
      box_1 = extruded_entities_1[i];
      break;
    }
  }
  const int number_of_layers = 3;
  auto extruded_entities_2 =
      extrude_entities({rectangle_2}, dr, true, number_of_layers);
  std::pair<int, int> box_2 = {-1, -1};
  for (std::size_t i = 0; i < extruded_entities_2.size(); i++) {
    if (extruded_entities_2[i].first == 3) {
      box_2 = extruded_entities_2[i];
      break;
    }
  }

  synchronize();

  set_granularity(0.5);
  generate_mesh(3, 1, 6);

  const auto nodes_1 = get_nodes_of<3>(box_1);
  bool at_least_one_node_not_on_plane = false;
  for (std::size_t i = 0; i < nodes_1.size(); i++) {
    if (fabs(double(int(nodes_1[i].position[2] /
                            (dr[2] / double(number_of_layers)) +
                        1e-10)) -
             (nodes_1[i].position[2] / (dr[2] / double(number_of_layers)))) >
        1e-10) {
      at_least_one_node_not_on_plane = true;
      break;
    }
  }
  ASSERT_THAT(at_least_one_node_not_on_plane, Eq(true));

  const auto nodes_2 = get_nodes_of<3>(box_2);
  for (std::size_t i = 0; i < nodes_2.size(); i++) {
    ASSERT_THAT(
        double(int(nodes_2[i].position[2] / (dr[2] / double(number_of_layers)) +
                   1e-10)),
        DoubleEq(nodes_2[i].position[2] / (dr[2] / double(number_of_layers))));
  }
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108