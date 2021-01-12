// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/tensor/Tensor.h"
#include "ae108/elements/tensor/as_matrix_of_columns.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace elements {
namespace tensor {
namespace {

struct as_matrix_of_columns_Test : Test {};

TEST_F(as_matrix_of_columns_Test, works_for_const_matrices) {
  const Tensor<int, 2, 3> x = {{{{7, 8, 9}}, {{10, 11, 12}}}};

  ASSERT_THAT(as_matrix_of_columns(&x).rows(), Eq(3));
  ASSERT_THAT(as_matrix_of_columns(&x).cols(), Eq(2));
  EXPECT_THAT(as_matrix_of_columns(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_matrix_of_columns(&x)(1, 0), Eq(8));
  EXPECT_THAT(as_matrix_of_columns(&x)(2, 0), Eq(9));
  EXPECT_THAT(as_matrix_of_columns(&x)(0, 1), Eq(10));
  EXPECT_THAT(as_matrix_of_columns(&x)(1, 1), Eq(11));
  EXPECT_THAT(as_matrix_of_columns(&x)(2, 1), Eq(12));
}

TEST_F(as_matrix_of_columns_Test, works_for_nonconst_matrices) {
  const Tensor<int, 2, 3> x = {{{{7, 8, 9}}, {{10, 11, 12}}}};

  ASSERT_THAT(as_matrix_of_columns(&x).rows(), Eq(3));
  ASSERT_THAT(as_matrix_of_columns(&x).cols(), Eq(2));
  EXPECT_THAT(as_matrix_of_columns(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_matrix_of_columns(&x)(1, 0), Eq(8));
  EXPECT_THAT(as_matrix_of_columns(&x)(2, 0), Eq(9));
  EXPECT_THAT(as_matrix_of_columns(&x)(0, 1), Eq(10));
  EXPECT_THAT(as_matrix_of_columns(&x)(1, 1), Eq(11));
  EXPECT_THAT(as_matrix_of_columns(&x)(2, 1), Eq(12));
}

TEST_F(as_matrix_of_columns_Test, modifying_nonconst_matrix_works) {
  Tensor<int, 2, 3> x = {{{{7, 8, 9}}, {{10, 11, 12}}}};

  as_matrix_of_columns (&x)(0, 0) = 3;

  EXPECT_THAT(x.at(0).at(0), Eq(3));
}

TEST_F(as_matrix_of_columns_Test, zero_columns_are_permitted) {
  const auto x = Tensor<int, 0, 2>();

  EXPECT_THAT(as_matrix_of_columns(&x).rows(), Eq(2));
  EXPECT_THAT(as_matrix_of_columns(&x).cols(), Eq(0));
}

TEST_F(as_matrix_of_columns_Test, one_row_is_permitted) {
  auto x = Tensor<int, 2, 1>();

  EXPECT_THAT(as_matrix_of_columns(&x).rows(), Eq(1));
  EXPECT_THAT(as_matrix_of_columns(&x).cols(), Eq(2));
}

TEST_F(as_matrix_of_columns_Test, one_column_is_permitted) {
  auto x = Tensor<int, 1, 2>();

  EXPECT_THAT(as_matrix_of_columns(&x).rows(), Eq(2));
  EXPECT_THAT(as_matrix_of_columns(&x).cols(), Eq(1));
}

TEST_F(as_matrix_of_columns_Test, one_row_is_permitted_in_const_case) {
  const auto x = Tensor<int, 2, 1>();

  EXPECT_THAT(as_matrix_of_columns(&x).rows(), Eq(1));
  EXPECT_THAT(as_matrix_of_columns(&x).cols(), Eq(2));
}

TEST_F(as_matrix_of_columns_Test, one_column_is_permitted_in_const_case) {
  const auto x = Tensor<int, 1, 2>();

  EXPECT_THAT(as_matrix_of_columns(&x).rows(), Eq(2));
  EXPECT_THAT(as_matrix_of_columns(&x).cols(), Eq(1));
}

} // namespace
} // namespace tensor
} // namespace elements
} // namespace ae108