// © 2021 ETH Zurich, Mechanics and Materials Lab
//
// This file is part of ae108.
//
// ae108 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any
// later version.
//
// ae108 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ae108. If not, see <https://www.gnu.org/licenses/>.

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