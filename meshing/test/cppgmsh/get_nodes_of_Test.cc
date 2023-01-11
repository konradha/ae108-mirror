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
#include "ae108/meshing/cppgmsh/fragment_entities.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Ne;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_nodes_Test : Test {};

TEST_F(get_nodes_Test, returns_correct_nodes) {
  const auto gmshContext = Context(0, 0);

  const auto box = construct_box({0., 0., 0.}, {1., 1., 1.});
  synchronize();

  set_granularity(1.);
  generate_mesh();

  const auto nodes = get_nodes_of<3>(box);
  ASSERT_THAT(nodes.size(), Eq(14));
  ASSERT_THAT(nodes[0].position,
              ElementsAre(DoubleEq(0), DoubleEq(0), DoubleEq(1)));
  ASSERT_THAT(nodes[1].position,
              ElementsAre(DoubleEq(0), DoubleEq(0), DoubleEq(0)));
  ASSERT_THAT(nodes[2].position,
              ElementsAre(DoubleEq(0), DoubleEq(1), DoubleEq(1)));
  ASSERT_THAT(nodes[3].position,
              ElementsAre(DoubleEq(0), DoubleEq(1), DoubleEq(0)));
  ASSERT_THAT(nodes[4].position,
              ElementsAre(DoubleEq(1), DoubleEq(0), DoubleEq(1)));
  ASSERT_THAT(nodes[5].position,
              ElementsAre(DoubleEq(1), DoubleEq(0), DoubleEq(0)));
  ASSERT_THAT(nodes[6].position,
              ElementsAre(DoubleEq(1), DoubleEq(1), DoubleEq(1)));
  ASSERT_THAT(nodes[7].position,
              ElementsAre(DoubleEq(1), DoubleEq(1), DoubleEq(0)));
}

TEST_F(get_nodes_Test, removes_duplicates) {
  const auto gmshContext = Context(0, 0);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({0., 0., 1.}, {1., 1., 1.});
  const auto fused_entities = fragment_entities({box_1, box_2});
  synchronize();

  set_granularity(1);
  generate_mesh(3, 1, 6);

  const auto duplicate_nodes = get_nodes_of<3>({3, -1}, false);
  const auto nodes = get_nodes_of<3>({3, -1}, true);

  ASSERT_THAT(duplicate_nodes.size(), Ne(nodes.size()));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108