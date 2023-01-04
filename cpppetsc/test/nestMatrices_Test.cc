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
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/createTransformInput.h"
#include "ae108/cpppetsc/multiply.h"
#include "ae108/cpppetsc/nestMatrices.h"
#include "ae108/cppptest/Matchers.h"
#include <cmath>
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::DoubleEq;
using testing::Pair;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct nestMatrices_Test : Test {
  using matrix_type = Matrix<Policy>;
  using vector_type = Vector<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(nestMatrices_Test, Policies);

TYPED_TEST(nestMatrices_Test, result_has_sum_of_rows) {
  using matrix_type = typename TestFixture::matrix_type;

  auto x = matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });
  auto y = matrix_type::fromList({
      {3., 3.},
  });

  const auto nested = nestMatrices<TypeParam>({
      {&x},
      {&y},
  });

  EXPECT_THAT(nested.size(), Pair(3, 2));
}

TYPED_TEST(nestMatrices_Test, result_has_sum_of_cols) {
  using matrix_type = typename TestFixture::matrix_type;

  auto x = matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });
  auto y = matrix_type::fromList({
      {3., 3.},
      {6., 6.},
  });

  const auto nested = nestMatrices<TypeParam>({
      {&x, &y},
  });

  EXPECT_THAT(nested.size(), Pair(2, 4));
}

TYPED_TEST(nestMatrices_Test, throws_if_inconsistent_rows) {
  using matrix_type = typename TestFixture::matrix_type;

  auto u = matrix_type::fromList({
      {1.},
  });
  auto v = matrix_type::fromList({
      {2.},
  });
  auto w = matrix_type::fromList({
      {3.},
  });

  EXPECT_THROW(nestMatrices<TypeParam>({
                   {&u, &v},
                   {&w},
               }),
               InvalidParametersException);
}

TYPED_TEST(nestMatrices_Test, result_has_correct_first_column) {
  using matrix_type = typename TestFixture::matrix_type;

  auto u = matrix_type::fromList({
      {1.},
  });
  auto v = matrix_type::fromList({
      {2.},
  });
  auto w = matrix_type::fromList({
      {3.},
  });
  auto x = matrix_type::fromList({
      {4.},
  });

  const auto nested = nestMatrices<TypeParam>({
      {&u, &v},
      {&w, &x},
  });

  auto in = createTransformInput(nested);
  in.unwrap().replace()(0) = 1.;
  const auto entry = multiply(nested, in);

  EXPECT_THAT(entry.unwrap(), ScalarEqIfLocal(0, 1.));
  EXPECT_THAT(entry.unwrap(), ScalarEqIfLocal(1, 3.));
}

TYPED_TEST(nestMatrices_Test, result_has_correct_second_column) {
  using matrix_type = typename TestFixture::matrix_type;

  auto u = matrix_type::fromList({
      {1.},
  });
  auto v = matrix_type::fromList({
      {2.},
  });
  auto w = matrix_type::fromList({
      {3.},
  });
  auto x = matrix_type::fromList({
      {4.},
  });

  const auto nested = nestMatrices<TypeParam>({
      {&u, &v},
      {&w, &x},
  });

  auto in = createTransformInput(nested);
  in.unwrap().replace()(1) = 1.;
  const auto entry = multiply(nested, in);

  EXPECT_THAT(entry.unwrap(), ScalarEqIfLocal(2, 1.));
  EXPECT_THAT(entry.unwrap(), ScalarEqIfLocal(4, 3.));
}

} // namespace cpppetsc
} // namespace ae108