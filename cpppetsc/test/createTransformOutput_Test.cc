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
#include "ae108/cpppetsc/createTransformOutput.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEq;
using testing::DoubleEq;
using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct createTransformOutput_Test : Test {
  using matrix_type = Matrix<Policy>;
  using size_type = typename Vector<Policy>::size_type;
  using value_type = typename Vector<Policy>::value_type;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(createTransformOutput_Test, Policies);

TYPED_TEST(createTransformOutput_Test, vector_has_correct_size) {
  using size_type = typename TestFixture::size_type;
  using matrix_type = typename TestFixture::matrix_type;

  const auto rows = size_type{3};
  const auto cols = size_type{7};
  const auto matrix = matrix_type(rows, cols);

  const auto result = createTransformOutput(matrix);
  EXPECT_THAT(result.unwrap(), SizeIs(rows));
}

TYPED_TEST(createTransformOutput_Test, vector_has_zero_norm) {
  using size_type = typename TestFixture::size_type;
  using matrix_type = typename TestFixture::matrix_type;

  const auto rows = size_type{3};
  const auto cols = size_type{7};
  const auto matrix = matrix_type(rows, cols);

  const auto result = createTransformOutput(matrix);
  EXPECT_THAT(result.unwrap().norm(), DoubleEq(0.));
}

TYPED_TEST(createTransformOutput_Test, multiplication_works) {
  using matrix_type = typename TestFixture::matrix_type;
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  const auto rows = size_type{3};
  const auto cols = size_type{7};
  auto matrix = matrix_type(rows, cols);
  matrix.setZero();

  auto in = createTransformInput(matrix);
  in.unwrap().setZero();
  auto out = createTransformOutput(matrix);
  out.unwrap().setZero();

  out.unwrap().addAx(matrix, in);
  EXPECT_THAT(out.unwrap().norm(), ScalarEq(value_type{0.}));
}

} // namespace cpppetsc
} // namespace ae108