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
TYPED_TEST(asTransposedMatrix_Test, first_column_is_nonhermitian_transpose) {
  using value_type = typename TestFixture::matrix_type::value_type;
  const auto matrix = asTransposedMatrix(&this->complex_matrix);

  auto input = createTransformOutput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, value_type{0, 1.}));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, value_type{0, 2.}));
}

TYPED_TEST(asTransposedMatrix_Test, second_column_is_nonhermitian_transpose) {
  using value_type = typename TestFixture::matrix_type::value_type;
  const auto matrix = asTransposedMatrix(&this->complex_matrix);

  auto input = createTransformOutput(matrix);
  input.unwrap().replace()(1) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, value_type{0, 3.}));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, value_type{0, 4.}));
}
#endif

} // namespace
} // namespace cpppetsc
} // namespace ae108