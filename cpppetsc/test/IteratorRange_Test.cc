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

#include "ae108/cpppetsc/IteratorRange.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace cpppetsc {

struct IteratorRange_Test : Test {
  std::vector<int> vector = {1, 2};
  IteratorRange<std::vector<int>::const_iterator> range{vector.begin(),
                                                        vector.end()};
};

TEST_F(IteratorRange_Test, begin_returns_begin) {
  EXPECT_THAT(range.begin(), Eq(vector.begin()));
}

TEST_F(IteratorRange_Test, end_returns_end) {
  EXPECT_THAT(range.end(), Eq(vector.end()));
}
} // namespace cpppetsc
} // namespace ae108
