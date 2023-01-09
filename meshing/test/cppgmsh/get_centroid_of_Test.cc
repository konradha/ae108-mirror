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
#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include "ae108/meshing/cppgmsh/get_centroid_of.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_centroid_of_Test : Test {};

TEST_F(get_centroid_of_Test, returns_correct_get_centroid_of_line) {
  const auto gmshContext = Context(0, 0);
  const auto centroid = get_centroid_of(
      {1, gmsh::model::occ::addLine(gmsh::model::occ::addPoint(0, 0, 0),
                                    gmsh::model::occ::addPoint(1, 0, 0))});
  ASSERT_THAT(centroid, ElementsAre(DoubleEq(0.5), DoubleEq(0), DoubleEq(0)));
}

TEST_F(get_centroid_of_Test, returns_correct_get_centroid_of_rectangle) {
  const auto gmshContext = Context(0, 0);
  const auto centroid = get_centroid_of(construct_rectangle({0, 0, 0}, {1, 1}));
  ASSERT_THAT(centroid, ElementsAre(DoubleEq(.5), DoubleEq(.5), DoubleEq(0)));
}

TEST_F(get_centroid_of_Test, returns_correct_get_centroid_of_box) {
  const auto gmshContext = Context(0, 0);
  const auto centroid = get_centroid_of(construct_box({0, 0, 0}, {1, 1, 1}));
  ASSERT_THAT(centroid, ElementsAre(DoubleEq(.5), DoubleEq(.5), DoubleEq(.5)));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108