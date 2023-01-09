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
#include "ae108/meshing/cppgmsh/construct_box.h"
#include "ae108/meshing/cppgmsh/get_entities_in.h"
#include "ae108/meshing/cppgmsh/set_physical_group_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_entities_in_Test : Test {};

TEST_F(get_entities_in_Test, returns_correct_vertices_in_box) {
  const auto gmshContext = Context(0, 0);
  gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1);
  gmsh::model::occ::synchronize();
  const auto points = get_entities_in({{0, 0, 0}, {1, 1, 1}}, 0);
  ASSERT_THAT(points[0].first, Eq(0));
  ASSERT_THAT(points.size(), Eq(8));
}

TEST_F(get_entities_in_Test, returns_correct_lines_in_box) {
  const auto gmshContext = Context(0, 0);
  gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1);
  gmsh::model::occ::synchronize();
  const auto lines = get_entities_in({{0, 0, 0}, {1, 1, 1}}, 1);
  ASSERT_THAT(lines[0].first, Eq(1));
  ASSERT_THAT(lines.size(), Eq(12));
}

TEST_F(get_entities_in_Test, returns_correct_surfaces_in_box) {
  const auto gmshContext = Context(0, 0);
  gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1);
  gmsh::model::occ::synchronize();
  const auto surfaces = get_entities_in({{0, 0, 0}, {1, 1, 1}}, 2);
  ASSERT_THAT(surfaces[0].first, Eq(2));
  ASSERT_THAT(surfaces.size(), Eq(6));
}

TEST_F(get_entities_in_Test, returns_correct_volumes_in_box) {
  const auto gmshContext = Context(0, 0);
  gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1);
  gmsh::model::occ::synchronize();
  const auto volumes = get_entities_in({{0, 0, 0}, {1, 1, 1}}, 3);
  ASSERT_THAT(volumes[0].first, Eq(3));
  ASSERT_THAT(volumes.size(), Eq(1));
}

TEST_F(get_entities_in_Test, returns_the_correct_entities) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({10., 10., 10.}, {1., 1., 1.});
  synchronize();

  const auto group = set_physical_group_of({box_1, box_2});
  const auto entities = get_entities_in(group);

  ASSERT_THAT(entities.size(), Eq(2));
  ASSERT_THAT(entities[0], Eq(box_1));
  ASSERT_THAT(entities[1], Eq(box_2));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108