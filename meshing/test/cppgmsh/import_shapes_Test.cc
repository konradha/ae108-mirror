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
#include "ae108/meshing/cppgmsh/fragment_entities.h"
#include "ae108/meshing/cppgmsh/import_shapes.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include "ae108/meshing/cppgmsh/write_file.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct import_shapes_Test : Test {};

TEST_F(import_shapes_Test, imports_box) {
  const auto gmshContext = Context(0, 0, true, false);

  const auto origin = std::array<double, 3>{0, 0, 0};
  const auto side_lengths = std::array<double, 3>{1, 2, 3};

  const auto box_0 = construct_box(origin, side_lengths);
  synchronize();
  write_file("box.stp");

  const auto imported_shapes = import_shapes("box.stp");
  ASSERT_THAT(imported_shapes.size(), Eq(1));
  const auto box_1 = imported_shapes[0];

  const auto fragments = fragment_entities({box_0, box_1});
  ASSERT_THAT(fragments.size(), Eq(1));
}

TEST_F(import_shapes_Test, imports_box_with_lower_dimensional_entities) {
  const auto gmshContext = Context(0, 0, true, false);

  const auto origin = std::array<double, 3>{0, 0, 0};
  const auto side_lengths = std::array<double, 3>{1, 2, 3};

  construct_box(origin, side_lengths);
  synchronize();
  write_file("box.stp");

  const auto imported_shapes = import_shapes("box.stp", false);
  std::array<bool, 4> includes_iD_entities = {false, false, false, false};
  for (std::size_t i = 0; i < imported_shapes.size(); i++) {
    for (int d = 0; d < 4; d++)
      if (imported_shapes[i].first == d)
        includes_iD_entities[d] = true;
  }
  for (std::size_t d = 0; d < 4; d++)
    ASSERT_THAT(includes_iD_entities[d], Eq(true));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108