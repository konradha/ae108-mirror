// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/copy.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct copy_Test : Test {
  using vector_type = Vector<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(copy_Test, Policies);

TYPED_TEST(copy_Test, copies_values) {
  using vector_type = typename TestFixture::vector_type;

  const auto from = tag<DistributedTag>(vector_type::fromList({3., 7.}));
  auto to = tag<DistributedTag>(vector_type(2));

  copy(from, &to);

  EXPECT_THAT(to.unwrap(), ScalarEqIfLocal(0, 3.));
  EXPECT_THAT(to.unwrap(), ScalarEqIfLocal(1, 7.));
}

} // namespace cpppetsc
} // namespace ae108