// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
