// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
//
// This file is part of ae108.
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

#include "ae108/assembly/utilities/resizeIfPossible.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <gmock/gmock.h>
#include <vector>

using testing::Eq;

namespace ae108 {
namespace assembly {
namespace utilities {
namespace {

TEST(resizeIfPossible_Test, resizes_vector) {
  std::vector<int> x;
  const std::vector<int>::size_type newSize = 77;
  resizeIfPossible(&x, newSize);
  EXPECT_THAT(x.size(), Eq(newSize));
}

TEST(resizeIfPossible_Test, does_not_resize_array) {
  constexpr std::size_t size = 2;
  auto x = std::array<int, 2>();
  resizeIfPossible(&x, size);
  EXPECT_THAT(x.size(), Eq(size));
}
} // namespace
} // namespace utilities
} // namespace assembly
} // namespace ae108
