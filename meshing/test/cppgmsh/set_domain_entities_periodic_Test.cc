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
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_periodic_nodes_of.h"
#include "ae108/meshing/cppgmsh/set_domain_entities_periodic.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Gt;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct set_domain_entities_periodic_Test : Test {};

TEST_F(set_domain_entities_periodic_Test, fuse_single_entity) {
  const auto gmshContext = Context(0, 0);
}

TEST_F(set_domain_entities_periodic_Test,
       returns_correct_mesh_of_periodic_rectangles) {
  const auto gmshContext = Context(0, 0);
  const auto source_rect = construct_rectangle({0, 0, 0}, {1, 1});
  const auto target_rect = construct_rectangle({0, 0, 1}, {1, 1});
  synchronize();

  const auto source_bbox =
      BoundingBox<std::array<double, 3>>{{0, 0, 0}, {1, 1, 0}};
  const auto target_bbox =
      BoundingBox<std::array<double, 3>>{{0, 0, 1}, {1, 1, 1}};

  set_domain_entities_periodic(source_bbox, target_bbox, 2);
  generate_mesh(2);

  std::pair<int, int> source;
  std::array<double, 16> transform;
  const auto [targetNodeTags, sourceNodeTags] =
      get_periodic_nodes_of(target_rect, &source, &transform);

  ASSERT_THAT(source, Eq(source_rect));
  ASSERT_THAT(targetNodeTags.size(), Gt(0));
  ASSERT_THAT(sourceNodeTags.size(), Eq(targetNodeTags.size()));
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