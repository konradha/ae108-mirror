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

template <class Policy> struct asInverseMatrix_Test : Test {
  using matrix_type = Matrix<Policy>;

  const matrix_type matrix = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(asInverseMatrix_Test, Policies);

TYPED_TEST(asInverseMatrix_Test, first_column_is_correct) {
  const auto matrix = asInverseMatrix(&this->matrix);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 1.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 0.));
}

TYPED_TEST(asInverseMatrix_Test, second_column_is_correct) {
  const auto matrix = asInverseMatrix(&this->matrix);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(1) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 0.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, .5));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108