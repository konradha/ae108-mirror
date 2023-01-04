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
#include "ae108/cpppetsc/asThAT.h"
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

template <class Policy> struct asThAT_Test : Test {
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
      {3.},
      {4.},
  });
#ifdef AE108_PETSC_COMPLEX
  const matrix_type complex_transform = matrix_type::fromList({
      {typename matrix_type::value_type{0, 3}, 0.},
      {0., typename matrix_type::value_type{0, 4}},
  });
#endif
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(asThAT_Test, Policies);

TYPED_TEST(asThAT_Test, first_column_is_correct_for_square_transform) {
  const auto matrix = asThAT(&this->matrix, &this->square_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 3. * 3. + 10. * 5.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 3. * 4. + 10. * 6.));
}

TYPED_TEST(asThAT_Test, second_column_is_correct_for_square_transform) {
  const auto matrix = asThAT(&this->matrix, &this->square_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(1) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 4. * 3. + 12. * 5.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 4. * 4. + 12. * 6.));
}

TYPED_TEST(asThAT_Test, result_is_correct_for_nonsquare_transform) {
  const auto matrix = asThAT(&this->matrix, &this->nonsquare_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 3. * 3. + 8. * 4.));
}

#ifdef AE108_PETSC_COMPLEX
TYPED_TEST(asThAT_Test, first_column_correct_for_complex_transform) {
  const auto matrix = asThAT(&this->matrix, &this->complex_transform);

  auto input = createTransformInput(matrix);
  input.unwrap().replace()(0) = 1.;
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 9.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 0.));
}

TYPED_TEST(asThAT_Test, second) {
  const auto matrix = asThAT(&this->matrix, &this->complex_transform);

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