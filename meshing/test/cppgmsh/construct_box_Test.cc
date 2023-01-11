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
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct construct_box_Test : Test {};

TEST_F(construct_box_Test, constructs_box) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto origin = std::array<double, 3>{0, 0, 0};
  const auto side_lengths = std::array<double, 3>{1, 2, 3};

  const auto box = construct_box(origin, side_lengths);

  ASSERT_THAT(box.first, Eq(3));
  ASSERT_THAT(box.second, Eq(1));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108