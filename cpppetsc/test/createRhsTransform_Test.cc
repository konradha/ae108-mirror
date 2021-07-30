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
#include "ae108/cpppetsc/createRhsTransform.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Eq;
using testing::Pair;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct createRhsTransform_Test : Test {
  using matrix_type = Matrix<Policy>;
  using size_type = typename matrix_type::size_type;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(createRhsTransform_Test, Policies);

TYPED_TEST(createRhsTransform_Test, has_correct_size) {
  using matrix_type = typename TestFixture::matrix_type;
  using size_type = typename TestFixture::size_type;

  const auto matrix = matrix_type(2, 3);

  const auto columns = size_type{7};
  const auto result = createRhsTransform(matrix, columns);

  EXPECT_THAT(result.size(), Pair(Eq(matrix.size().second), Eq(columns)));
}

TYPED_TEST(createRhsTransform_Test, can_be_multiplied) {
  using matrix_type = typename TestFixture::matrix_type;
  using size_type = typename TestFixture::size_type;

  auto matrix = matrix_type(2, 3);
  matrix.setZero();

  const auto columns = size_type{7};
  auto transform = createRhsTransform(matrix, columns);
  transform.setZero();

  const auto result = matrix_type::fromProduct(matrix, transform);
  EXPECT_THAT(result.norm(), DoubleEq(0.));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108