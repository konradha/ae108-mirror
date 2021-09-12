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
#include "ae108/cpppetsc/asTransposedMatrix.h"
#include "ae108/cpppetsc/createTransformOutput.h"
#include "ae108/cpppetsc/multiply.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct asTransposedMatrix_Test : Test {
  using matrix_type = Matrix<Policy>;

  const matrix_type matrix = matrix_type::fromList({
      {1., 2.},
      {3., 4.},
  });
#ifdef AE108_PETSC_COMPLEX
  const matrix_type complex_matrix = matrix_type::fromList({
      {
          typename matrix_type::value_type{0, 1.},
          typename matrix_type::value_type{0, 2.},
      },
      {
          typename matrix_type::value_type{0, 3.},
          typename matrix_type::value_type{0, 4.},
      },
  });
#endif
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(asTransposedMatrix_Test, Policies);

TYPED_TEST(asTransposedMatrix_Test, first_column_is_correct) {
  const auto matrix = asTransposedMatrix(&this->matrix);

  auto input = createTransformOutput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 1.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 2.));
}

TYPED_TEST(asTransposedMatrix_Test, second_column_is_correct) {
  const auto matrix = asTransposedMatrix(&this->matrix);

  auto input = createTransformOutput(matrix);
  input.unwrap().replace()(1) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 3.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 4.));
}

#ifdef AE108_PETSC_COMPLEX
TYPED_TEST(asTransposedMatrix_Test, first_column_is_hermitian_transpose) {
  using value_type = typename TestFixture::matrix_type::value_type;
  const auto matrix = asTransposedMatrix(&this->complex_matrix);

  auto input = createTransformOutput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, value_type{0, -1.}));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, value_type{0, -2.}));
}

TYPED_TEST(asTransposedMatrix_Test, second_column_is_hermitian_transpose) {
  using value_type = typename TestFixture::matrix_type::value_type;
  const auto matrix = asTransposedMatrix(&this->complex_matrix);

  auto input = createTransformOutput(matrix);
  input.unwrap().replace()(1) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, value_type{0, -3.}));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, value_type{0, -4.}));
}
#endif

} // namespace
} // namespace cpppetsc
} // namespace ae108