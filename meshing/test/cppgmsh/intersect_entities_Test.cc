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
#include "ae108/meshing/cppgmsh/get_entities_of.h"
#include "ae108/meshing/cppgmsh/intersect_entities.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct intersect_entities_Test : Test {};

TEST_F(intersect_entities_Test, self_intersect_single_entitity) {
  const auto gmshContext = Context(0, 0);

  const auto object = construct_box({0, 0, 0}, {1, 1, 1});

  const auto result = intersect_entities({object}, {object});
  ASSERT_THAT(result.size(), Eq(1));
}

TEST_F(intersect_entities_Test, intersect_two_entitities) {
  const auto gmshContext = Context(0, 0);

  const auto object = construct_box({0, 0, 0}, {1, 1, 1});
  const auto tool = construct_box({.5, .5, .5}, {1, 1, 1});

  const auto result = intersect_entities({object}, {tool});
  ASSERT_THAT(result.size(), Eq(1));
}

TEST_F(intersect_entities_Test, intersect_unconnected_entitities) {
  const auto gmshContext = Context(0, 0);

  const auto object = construct_box({0, 0, 0}, {1, 1, 1});
  const auto tool = construct_box({2, 2, 2}, {1, 1, 1});

  const auto result = intersect_entities({object}, {tool});
  ASSERT_THAT(result.size(), Eq(0));
}

TEST_F(intersect_entities_Test, keep_object) {
  const auto gmshContext = Context(0, 0);

  const auto object = construct_box({0, 0, 0}, {1, 1, 1});
  const auto tool = construct_box({.5, .5, .5}, {1, 1, 1});

  intersect_entities({object}, {tool}, false, true);
  synchronize();
  const auto entities = get_entities_of(3);
  ASSERT_THAT(entities.size(), Eq(2));
}

TEST_F(intersect_entities_Test, keep_tool) {
  const auto gmshContext = Context(0, 0);

  const auto object = construct_box({0, 0, 0}, {1, 1, 1});
  const auto tool = construct_box({.5, .5, .5}, {1, 1, 1});

  intersect_entities({object}, {tool}, true, false);
  synchronize();
  const auto entities = get_entities_of(3);
  ASSERT_THAT(entities.size(), Eq(2));
}

TEST_F(intersect_entities_Test, keep_object_and_tool) {
  const auto gmshContext = Context(0, 0);

  const auto object = construct_box({0, 0, 0}, {1, 1, 1});
  const auto tool = construct_box({.5, .5, .5}, {1, 1, 1});

  intersect_entities({object}, {tool}, false, false);
  synchronize();
  const auto entities = get_entities_of(3);
  ASSERT_THAT(entities.size(), Eq(3));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108