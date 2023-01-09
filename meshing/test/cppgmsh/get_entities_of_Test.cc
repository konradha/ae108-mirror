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
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_entities_of_Test : Test {};

TEST_F(get_entities_of_Test, gets_enities_of_dim3) {
  const auto gmshContext = Context(0, 0, true, true);

  std::pair<int, int> box = construct_box(std::array<double, 3>{0., 0., 0.},
                                          std::array<double, 3>{1., 1., 1.});
  synchronize();

  const auto entities = get_entities_of(3);
  ASSERT_THAT(entities.size(), Eq(1));
  ASSERT_THAT(entities[0].first, Eq(3));
  ASSERT_THAT(entities[0].second, Eq(box.second));
}

TEST_F(get_entities_of_Test, gets_enities_of_dim2) {
  const auto gmshContext = Context(0, 0, true, true);

  construct_box(std::array<double, 3>{0., 0., 0.},
                std::array<double, 3>{1., 1., 1.});
  synchronize();

  const auto entities = get_entities_of(2);
  ASSERT_THAT(entities.size(), Eq(6));
  ASSERT_THAT(entities[0].first, Eq(2));
}

TEST_F(get_entities_of_Test, gets_enities_of_all_dimensions) {
  const auto gmshContext = Context(0, 0, true, true);

  construct_box(std::array<double, 3>{0., 0., 0.},
                std::array<double, 3>{1., 1., 1.});
  synchronize();

  const auto entities = get_entities_of();
  ASSERT_THAT(entities.size(),
              Eq(1 /*volume*/ + 6 /*surfaces*/ + 12 /*edges*/ + 8 /*points*/));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108