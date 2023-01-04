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

template <class Policy> struct multiply_Test : Test {
  using matrix_type = Matrix<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(multiply_Test, Policies);

TYPED_TEST(multiply_Test, multiplies_matrix_with_vector_correctly) {
  const auto matrix = TestFixture::matrix_type::fromList({
      {1., 2.},
      {3., 4.},
  });

  auto input = createTransformInput(matrix);
  {
    auto replacer = input.unwrap().replace();
    replacer(0) = 1.;
    replacer(1) = -2.;
  }
  const auto result = multiply(matrix, input);

  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, -3.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, -5.));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108