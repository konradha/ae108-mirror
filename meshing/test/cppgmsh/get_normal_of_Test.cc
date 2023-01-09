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
#include "ae108/meshing/cppgmsh/get_normal_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_normal_of_Test : Test {};

TEST_F(get_normal_of_Test, returns_correct_normal_of_rectangle) {
  const auto gmshContext = Context(0, 0);
  const auto rectangle = construct_rectangle({0, 0, 0}, {1, 1});
  synchronize();
  const auto normal = get_normal_of(rectangle.second);
  ASSERT_THAT(normal, ElementsAre(DoubleEq(0), DoubleEq(0), DoubleEq(1)));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108