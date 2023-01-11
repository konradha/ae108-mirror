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

#include "ae108/meshing/cppgmsh/as_affine_transform.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct as_affine_transform_Test : Test {};

TEST_F(as_affine_transform_Test, returns_correct_translation_transform) {
  const auto transform = as_affine_transform({1., 2., 3.});
  ASSERT_THAT(transform, ElementsAre(DoubleEq(1.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(1.), DoubleEq(0.), DoubleEq(1.),
                                     DoubleEq(0.), DoubleEq(2.), DoubleEq(0.),
                                     DoubleEq(0.), DoubleEq(1.), DoubleEq(3.),
                                     DoubleEq(0.), DoubleEq(0.), DoubleEq(0.),
                                     DoubleEq(1.)));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108