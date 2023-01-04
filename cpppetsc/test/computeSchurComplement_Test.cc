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
#include "ae108/cpppetsc/asSchurComplement.h"
#include "ae108/cpppetsc/computeSchurComplement.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct computeSchurComplement_Test : Test {
  using matrix_type = Matrix<Policy>;

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

  const matrix_type mat = matrix_type::fromList({
      {1., 0., 3., 0.},
      {0., 2., 0., 4.},
      {5., 0., 7., 0.},
      {0., 6., 0., 8.},
  });
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(computeSchurComplement_Test, Policies);

TYPED_TEST(computeSchurComplement_Test,
           schur_complement_created_using_submatrices_has_correct_elements) {
  const auto matrix = computeSchurComplement(asSchurComplement(
      &this->mat_00, &this->mat_01, &this->mat_10, &this->mat_11));

  EXPECT_THAT(matrix, ScalarEqIfLocal(0, 0, 7. - 3. * 5. / 1.));
  EXPECT_THAT(matrix, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(matrix, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(matrix, ScalarEqIfLocal(1, 1, 8. - 4. * 6. / 2.));
}

TYPED_TEST(computeSchurComplement_Test,
           schur_complement_created_using_indices_has_correct_elements) {
  const auto matrix =
      computeSchurComplement(asSchurComplement(&this->mat, {2, 3}));

  EXPECT_THAT(matrix, ScalarEqIfLocal(0, 0, 7. - 3. * 5. / 1.));
  EXPECT_THAT(matrix, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(matrix, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(matrix, ScalarEqIfLocal(1, 1, 8. - 4. * 6. / 2.));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108