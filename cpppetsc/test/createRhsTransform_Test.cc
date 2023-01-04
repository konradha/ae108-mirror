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

#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/createRhsTransform.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Eq;
using testing::Pair;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct createRhsTransform_Test : Test {
  using matrix_type = Matrix<Policy>;
  using size_type = typename matrix_type::size_type;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(createRhsTransform_Test, Policies);

TYPED_TEST(createRhsTransform_Test, has_correct_global_size) {
  using matrix_type = typename TestFixture::matrix_type;
  using size_type = typename TestFixture::size_type;

  const auto matrix = matrix_type(2, 3);

  const auto columns = size_type{7};
  const auto result = createRhsTransform(matrix, columns);

  EXPECT_THAT(result.size(), Pair(Eq(matrix.size().second), Eq(columns)));
}

TYPED_TEST(createRhsTransform_Test, has_correct_local_size) {
  using matrix_type = typename TestFixture::matrix_type;
  using size_type = typename TestFixture::size_type;

  const auto matrix = matrix_type(2, 3);

  const auto columns = size_type{7};
  const auto result =
      createRhsTransform(matrix, typename matrix_type::LocalCols{columns},
                         typename matrix_type::GlobalCols{PETSC_DETERMINE});

  EXPECT_THAT(result.localSize(),
              Pair(Eq(matrix.localSize().second), Eq(columns)));
}

TYPED_TEST(createRhsTransform_Test, can_be_multiplied) {
  using matrix_type = typename TestFixture::matrix_type;
  using size_type = typename TestFixture::size_type;

  auto matrix = matrix_type(2, 3);
  matrix.setZero();

  const auto columns = size_type{7};
  auto transform = createRhsTransform(matrix, columns);
  transform.setZero();

  const auto result = matrix_type::fromProduct(matrix, transform);
  EXPECT_THAT(result.norm(), DoubleEq(0.));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108