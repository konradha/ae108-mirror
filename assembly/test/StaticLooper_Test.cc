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

#include "ae108/assembly/utilities/StaticLooper.h"
#include <cstddef>
#include <gmock/gmock.h>
#include <vector>

using testing::ElementsAre;
using testing::IsEmpty;
using testing::Test;

namespace ae108 {
namespace assembly {
namespace utilities {
namespace {

struct Group {
  template <std::size_t N> std::size_t get() const { return N; }
};

struct StaticLooper_Test : public Test {};

TEST_F(StaticLooper_Test, check_loop_to_two) {
  std::vector<std::size_t> result;
  Group group;
  StaticLooper<2>::run(
      [&result](const std::size_t value) { result.push_back(value); }, group);

  EXPECT_THAT(result, ElementsAre(0, 1));
}

TEST_F(StaticLooper_Test, check_empty_loop) {
  std::vector<std::size_t> result;
  Group group;
  StaticLooper<0>::run(
      [&result](const std::size_t value) { result.push_back(value); }, group);

  EXPECT_THAT(result, IsEmpty());
}

TEST_F(StaticLooper_Test, loop_with_parameter) {
  std::vector<std::size_t> result;
  Group group;
  StaticLooper<2>::run(
      [&result](const std::size_t value, const std::size_t parameter) {
        result.push_back(value + parameter);
      },
      group, 7);

  EXPECT_THAT(result, ElementsAre(7, 8));
}
} // namespace
} // namespace utilities
} // namespace assembly
} // namespace ae108
