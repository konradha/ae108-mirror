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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/clone.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::AlmostEqIfLocal;
using testing::Pair;
using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct clone_Test : Test {
  using vector_type = Vector<Policy>;
  using matrix_type = Matrix<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(clone_Test, Policies);

TYPED_TEST(clone_Test, clones_vectors) {
  using vector_type = typename TestFixture::vector_type;

  const auto from = tag<DistributedTag>(vector_type::fromList({3., 7.}));

  const auto result = clone(from);

  ASSERT_THAT(result.unwrap(), SizeIs(2));
  EXPECT_THAT(result.unwrap(), AlmostEqIfLocal(0, 3.));
  EXPECT_THAT(result.unwrap(), AlmostEqIfLocal(1, 7.));
}

TYPED_TEST(clone_Test, clones_matrices) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto from = matrix_type::fromList({{1., 2.}, {3., 4.}});

  const auto result = clone(from);

  ASSERT_THAT(result.size(), Pair(2, 2));
  EXPECT_THAT(result, AlmostEqIfLocal(0, 0, 1.));
  EXPECT_THAT(result, AlmostEqIfLocal(0, 1, 2.));
  EXPECT_THAT(result, AlmostEqIfLocal(1, 0, 3.));
  EXPECT_THAT(result, AlmostEqIfLocal(1, 1, 4.));
}

} // namespace cpppetsc
} // namespace ae108