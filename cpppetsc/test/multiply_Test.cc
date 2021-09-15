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
#include "ae108/cpppetsc/createTransformInput.h"
#include "ae108/cpppetsc/multiply.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct multiply_Test : Test {
  using matrix_type = Matrix<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(multiply_Test, Policies);

TYPED_TEST(multiply_Test, multiplies_matrix_with_vector_correctly) {
  const auto matrix = TestFixture::matrix_type::fromList({
      {1., 2.},
      {3., 4.},
  });

  auto input = createTransformInput(matrix);
  {
    auto replacer = input.unwrap().replace();
    replacer(0) = 1.;
    replacer(1) = -2.;
  }
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, -3.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, -5.));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108