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
#include "ae108/cpppetsc/asInverseMatrix.h"
#include "ae108/cpppetsc/asSchurComplement.h"
#include "ae108/cpppetsc/asThAT.h"
#include "ae108/cpppetsc/asTransposedMatrix.h"
#include "ae108/cpppetsc/computeElementsOfMatrix.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct computeElementsOfMatrix_Test : Test {
  using matrix_type = Matrix<Policy>;

  const matrix_type matrix = matrix_type::fromList({
      {1., 2.},
      {3., 4.},
  });

  const matrix_type transform = matrix_type::fromList({
      {-1., 2.},
      {0., -3.},
  });
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(computeElementsOfMatrix_Test, Policies);

TYPED_TEST(computeElementsOfMatrix_Test,
           computing_the_elements_of_transposed_matrix_works) {
  const auto result =
      computeElementsOfMatrix(asTransposedMatrix(&this->matrix));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 1.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 3.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 0, 2.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 1, 4.));
}

TYPED_TEST(computeElementsOfMatrix_Test,
           computing_the_elements_of_transformed_matrix_works) {
  const auto result =
      computeElementsOfMatrix(asThAT(&this->matrix, &this->transform));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 1.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 4.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 0, 7.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 1, 10.));
}

TYPED_TEST(computeElementsOfMatrix_Test,
           computing_the_elements_of_inverse_matrix_works) {
  const auto result = computeElementsOfMatrix(asInverseMatrix(&this->matrix));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, -2.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 1.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 0, 1.5));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 1, -.5));
}

TYPED_TEST(computeElementsOfMatrix_Test,
           computing_the_elements_of_default_matrix_works) {
  const auto result = computeElementsOfMatrix(this->matrix);

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 1.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 2.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 0, 3.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 1, 4.));
}

TYPED_TEST(computeElementsOfMatrix_Test,
           computing_the_elements_of_schur_complement_works) {
  using matrix_type = typename TestFixture::matrix_type;

  const matrix_type mat_00 = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  const matrix_type mat_01 = matrix_type::fromList({
      {3., 0.},
      {0., 4.},
  });
  const matrix_type mat_10 = matrix_type::fromList({
      {5., 0.},
      {0., 6.},
  });
  const matrix_type mat_11 = matrix_type::fromList({
      {7., 0.},
      {0., 8.},
  });

  const auto result = computeElementsOfMatrix(
      asSchurComplement(&mat_00, &mat_01, &mat_10, &mat_11));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, -8.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 1, -4.));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108