// © 2022 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/meshing/cppgmsh/remove_entities.h"
#include "ae108/meshing/cppgmsh/synchronize.h"

#include <gmock/gmock.h>

using testing::Eq;
using testing::Gt;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct remove_entities_Test : Test {};

TEST_F(remove_entities_Test, remove_box_keep_boundaries) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box = construct_box({0., 0., 0.}, {1., 1., 1.});

  remove_entities({box});

  synchronize();
  auto entities = get_entities_of(3);
  ASSERT_THAT(entities.size(), Eq(0));
  entities = get_entities_of();
  ASSERT_THAT(entities.size(), Gt(0));
}

TEST_F(remove_entities_Test, remove_box_and_boundaries) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto box = construct_box({0., 0., 0.}, {1., 1., 1.});

  remove_entities({box}, true);

  const auto entities = get_entities_of();
  ASSERT_THAT(entities.size(), Eq(0));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108