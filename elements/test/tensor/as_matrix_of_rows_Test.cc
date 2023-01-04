// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/tensor/Tensor.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace elements {
namespace tensor {
namespace {

struct as_matrix_of_rows_Test : Test {};

TEST_F(as_matrix_of_rows_Test, works_for_const_matrices) {
  const Tensor<int, 2, 3> x = {{{{7, 8, 9}}, {{10, 11, 12}}}};

  ASSERT_THAT(as_matrix_of_rows(&x).rows(), Eq(2));
  ASSERT_THAT(as_matrix_of_rows(&x).cols(), Eq(3));
  EXPECT_THAT(as_matrix_of_rows(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_matrix_of_rows(&x)(0, 1), Eq(8));
  EXPECT_THAT(as_matrix_of_rows(&x)(0, 2), Eq(9));
  EXPECT_THAT(as_matrix_of_rows(&x)(1, 0), Eq(10));
  EXPECT_THAT(as_matrix_of_rows(&x)(1, 1), Eq(11));
  EXPECT_THAT(as_matrix_of_rows(&x)(1, 2), Eq(12));
}

TEST_F(as_matrix_of_rows_Test, works_for_nonconst_matrices) {
  const Tensor<int, 2, 3> x = {{{{7, 8, 9}}, {{10, 11, 12}}}};

  ASSERT_THAT(as_matrix_of_rows(&x).rows(), Eq(2));
  ASSERT_THAT(as_matrix_of_rows(&x).cols(), Eq(3));
  EXPECT_THAT(as_matrix_of_rows(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_matrix_of_rows(&x)(0, 1), Eq(8));
  EXPECT_THAT(as_matrix_of_rows(&x)(0, 2), Eq(9));
  EXPECT_THAT(as_matrix_of_rows(&x)(1, 0), Eq(10));
  EXPECT_THAT(as_matrix_of_rows(&x)(1, 1), Eq(11));
  EXPECT_THAT(as_matrix_of_rows(&x)(1, 2), Eq(12));
}

TEST_F(as_matrix_of_rows_Test, modifying_nonconst_matrix_works) {
  Tensor<int, 2, 3> x = {{{{7, 8, 9}}, {{10, 11, 12}}}};

  as_matrix_of_rows (&x)(0, 0) = 3;

  EXPECT_THAT(x.at(0).at(0), Eq(3));
}

TEST_F(as_matrix_of_rows_Test, zero_rows_are_permitted) {
  const auto x = Tensor<int, 0, 2>();

  EXPECT_THAT(as_matrix_of_rows(&x).rows(), Eq(0));
  EXPECT_THAT(as_matrix_of_rows(&x).cols(), Eq(2));
}

TEST_F(as_matrix_of_rows_Test, one_row_is_permitted) {
  auto x = Tensor<int, 1, 2>();

  EXPECT_THAT(as_matrix_of_rows(&x).rows(), Eq(1));
  EXPECT_THAT(as_matrix_of_rows(&x).cols(), Eq(2));
}

TEST_F(as_matrix_of_rows_Test, one_column_is_permitted) {
  auto x = Tensor<int, 2, 1>();

  EXPECT_THAT(as_matrix_of_rows(&x).rows(), Eq(2));
  EXPECT_THAT(as_matrix_of_rows(&x).cols(), Eq(1));
}

TEST_F(as_matrix_of_rows_Test, one_row_is_permitted_in_const_case) {
  const auto x = Tensor<int, 1, 2>();

  EXPECT_THAT(as_matrix_of_rows(&x).rows(), Eq(1));
  EXPECT_THAT(as_matrix_of_rows(&x).cols(), Eq(2));
}

TEST_F(as_matrix_of_rows_Test, one_column_is_permitted_in_const_case) {
  const auto x = Tensor<int, 2, 1>();

  EXPECT_THAT(as_matrix_of_rows(&x).rows(), Eq(2));
  EXPECT_THAT(as_matrix_of_rows(&x).cols(), Eq(1));
}

} // namespace
} // namespace tensor
} // namespace elements
} // namespace ae108