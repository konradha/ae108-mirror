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
#include "ae108/cpppetsc/asTransformedMatrix.h"
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

template <class Policy> struct asTransformedMatrix_Test : Test {
  using matrix_type = Matrix<Policy>;

  const matrix_type matrix = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  const matrix_type square_transform = matrix_type::fromList({
      {3., 4.},
      {5., 6.},
  });
  const matrix_type nonsquare_transform = matrix_type::fromList({
      {3., 4.},
  });
#ifdef AE108_PETSC_COMPLEX
  const matrix_type complex_transform = matrix_type::fromList({
      {typename matrix_type::value_type{0, 3}, 0.},
      {0., typename matrix_type::value_type{0, 4}},
  });
#endif
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(asTransformedMatrix_Test, Policies);

TYPED_TEST(asTransformedMatrix_Test,
           first_column_is_correct_for_square_transform) {
  const auto matrix =
      asTransformedMatrix(&this->matrix, &this->square_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 3. * 3. + 8. * 4.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 3. * 5. + 8. * 6.));
}

TYPED_TEST(asTransformedMatrix_Test,
           second_column_is_correct_for_square_transform) {
  const auto matrix =
      asTransformedMatrix(&this->matrix, &this->square_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(1) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 5. * 3. + 12. * 4.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 5. * 5. + 12. * 6.));
}

TYPED_TEST(asTransformedMatrix_Test,
           result_is_correct_for_nonsquare_transform) {
  const auto matrix =
      asTransformedMatrix(&this->matrix, &this->nonsquare_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 3. * 3. + 8. * 4.));
}

#ifdef AE108_PETSC_COMPLEX
TYPED_TEST(asTransformedMatrix_Test,
           first_column_correct_for_complex_transform) {
  const auto matrix =
      asTransformedMatrix(&this->matrix, &this->complex_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 9.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 0.));
}

TYPED_TEST(asTransformedMatrix_Test, second) {
  const auto matrix =
      asTransformedMatrix(&this->matrix, &this->complex_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(1) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 0.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 32.));
}
#endif

} // namespace
} // namespace cpppetsc
} // namespace ae108