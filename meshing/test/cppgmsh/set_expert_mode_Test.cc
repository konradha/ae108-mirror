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
#include "ae108/meshing/cppgmsh/set_expert_mode.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::DoubleEq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct set_expert_mode_Test : Test {};

TEST_F(set_expert_mode_Test, enables_expert_mode) {
  const auto gmshContext = Context(0, 0, true, true);

  set_expert_mode();
  double value;
  gmsh::option::getNumber("General.ExpertMode", value);
  ASSERT_THAT(value, DoubleEq(1.0));
}

TEST_F(set_expert_mode_Test, disables_expert_mode) {
  const auto gmshContext = Context(0, 0, true, true);

  set_expert_mode(false);
  double value;
  gmsh::option::getNumber("General.ExpertMode", value);
  ASSERT_THAT(value, DoubleEq(0.0));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108