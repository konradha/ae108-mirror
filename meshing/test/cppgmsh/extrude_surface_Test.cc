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
#include "ae108/meshing/cppgmsh/extrude_surface.h"
#include "ae108/meshing/cppgmsh/get_entities_of.h"
#include "ae108/meshing/cppgmsh/remove_entities.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include "ae108/meshing/cppgmsh/write_file.h"

#include <cmath>
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::DoubleEq;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct extrude_Test : Test {};

TEST_F(extrude_Test, extrude_rectangle_along_axis) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto rectangle = construct_rectangle({-.5, -.5, 0}, {1., 1.});
  const std::array<double, 3> from = {0, 0, 0};
  const std::array<double, 3> to = {0, 0, 10};

  extrude_surface(rectangle.second, from, to);
  remove_entities({rectangle}, true);
  synchronize();

  const auto points = get_entities_of(0);
  ASSERT_THAT(points.size(), Eq(8));
  const auto edges = get_entities_of(1);
  ASSERT_THAT(edges.size(), Eq(12));
  const auto surfaces = get_entities_of(2);
  ASSERT_THAT(surfaces.size(), Eq(6));
  const auto volumes = get_entities_of(3);
  ASSERT_THAT(volumes.size(), Eq(1));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108