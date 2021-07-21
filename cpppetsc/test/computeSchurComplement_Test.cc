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
#include "ae108/cpppetsc/asSchurComplement.h"
#include "ae108/cpppetsc/computeSchurComplement.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::AlmostEqIfLocal;
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

  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 0, 7. - 3. * 5. / 1.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 0, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 1, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 1, 8. - 4. * 6. / 2.));
}

TYPED_TEST(computeSchurComplement_Test,
           schur_complement_created_using_indices_has_correct_elements) {
  const auto matrix =
      computeSchurComplement(asSchurComplement(&this->mat, {2, 3}));

  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 0, 7. - 3. * 5. / 1.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 0, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 1, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 1, 8. - 4. * 6. / 2.));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108