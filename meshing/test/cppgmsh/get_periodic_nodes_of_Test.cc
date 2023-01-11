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
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_periodic_nodes_of.h"
#include "ae108/meshing/cppgmsh/set_periodic.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Gt;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_periodic_nodes_of_Test : Test {};

TEST_F(get_periodic_nodes_of_Test, returns_correct_mesh_of_periodic_lines) {
  const auto gmshContext = Context(0, 0);
  const auto source_line = gmsh::model::occ::addLine(
      gmsh::model::occ::addPoint(0, 0, 0), gmsh::model::occ::addPoint(1, 0, 0));
  const auto target_line = gmsh::model::occ::addLine(
      gmsh::model::occ::addPoint(0, 0, 1), gmsh::model::occ::addPoint(1, 0, 1));
  gmsh::model::occ::synchronize();

  set_periodic({{1, target_line}}, {1, source_line},
               {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1});
  generate_mesh(1);

  std::array<double, 16> transform;
  std::pair<int, int> source_entity;
  const auto node_tags =
      get_periodic_nodes_of({1, target_line}, &source_entity, &transform);

  ASSERT_THAT(source_entity.first, Eq(1));
  ASSERT_THAT(source_entity.second, Eq(source_line));
  ASSERT_THAT(node_tags.first.size(), Gt(0));
  ASSERT_THAT(node_tags.second.size(), Gt(0));
  ASSERT_THAT(node_tags.second.size(), Eq(node_tags.first.size()));
  ASSERT_THAT(transform, ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(0.), DoubleEq(0.), DoubleEq(1.),
                                     DoubleEq(0.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(0.), DoubleEq(1.), DoubleEq(1.),
                                     DoubleEq(0.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(1.)));
}

TEST_F(get_periodic_nodes_of_Test,
       returns_correct_mesh_of_periodic_rectangles) {
  const auto gmshContext = Context(0, 0);
  const auto source_rect = gmsh::model::occ::addRectangle(0, 0, 0, 1, 1);
  const auto target_rect = gmsh::model::occ::addRectangle(0, 0, 1, 1, 1);
  gmsh::model::occ::synchronize();

  set_periodic({{2, target_rect}}, {2, source_rect},
               {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1});
  generate_mesh(2);

  std::array<double, 16> transform;
  std::pair<int, int> source_entity;
  const auto node_tags =
      get_periodic_nodes_of({2, target_rect}, &source_entity, &transform);

  ASSERT_THAT(source_entity.first, Eq(2));
  ASSERT_THAT(source_entity.second, Eq(source_rect));
  ASSERT_THAT(node_tags.first.size(), Gt(0));
  ASSERT_THAT(node_tags.second.size(), Gt(0));
  ASSERT_THAT(node_tags.second.size(), Eq(node_tags.first.size()));
  ASSERT_THAT(transform, ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(0.), DoubleEq(0.), DoubleEq(1.),
                                     DoubleEq(0.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(0.), DoubleEq(1.), DoubleEq(1.),
                                     DoubleEq(0.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(1.)));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108