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
#include "ae108/meshing/cppgmsh/copy_entities.h"
#include "ae108/meshing/cppgmsh/get_coords_of.h"
#include "ae108/meshing/cppgmsh/get_entities_of.h"
#include "ae108/meshing/cppgmsh/get_points_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"

#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct copy_entities_Test : Test {};

TEST_F(copy_entities_Test, copies_box) {
  const auto gmshContext = Context(0, 0, true, true);

  auto box = construct_box(std::array<double, 3>{0., 0., 0.},
                           std::array<double, 3>{1., 1., 1.});

  copy_entities({box});

  synchronize();
  const auto entities = get_entities_of(3);

  ASSERT_THAT(entities.size(), Eq(2));
  auto points_0 = get_points_of(entities[0]);
  auto points_1 = get_points_of(entities[1]);
  ASSERT_THAT(points_0.size(), Eq(points_1.size()));
  for (std::size_t i = 0; i < points_0.size(); i++) {
    auto coords_0 = get_coords_of(points_0[i]);
    auto coords_1 = get_coords_of(points_1[i]);
    for (std::size_t d = 0; d < 3; d++) {
      ASSERT_THAT(coords_0[d], DoubleEq(coords_1[d]));
    }
  }
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108