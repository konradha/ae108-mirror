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
#include "ae108/elements/tensor/as_vector.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace elements {
namespace tensor {
namespace {

struct as_vector_Test : Test {};

TEST_F(as_vector_Test, works_for_const_vectors) {
  const tensor::Tensor<int, 3> x = {{7, 8, 9}};

  ASSERT_THAT(as_vector(&x).rows(), Eq(3));
  ASSERT_THAT(as_vector(&x).cols(), Eq(1));
  EXPECT_THAT(as_vector(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_vector(&x)(1, 0), Eq(8));
  EXPECT_THAT(as_vector(&x)(2, 0), Eq(9));
}

TEST_F(as_vector_Test, works_for_nonconst_vectors) {
  tensor::Tensor<int, 3> x = {{7, 8, 9}};

  ASSERT_THAT(as_vector(&x).rows(), Eq(3));
  ASSERT_THAT(as_vector(&x).cols(), Eq(1));
  EXPECT_THAT(as_vector(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_vector(&x)(1, 0), Eq(8));
  EXPECT_THAT(as_vector(&x)(2, 0), Eq(9));
}

TEST_F(as_vector_Test, modifying_nonconst_vector_works) {
  tensor::Tensor<int, 3> x = {{7, 8, 9}};

  as_vector (&x)(0, 0) = 3;

  EXPECT_THAT(x.at(0), Eq(3));
}
TEST_F(as_vector_Test, works_for_const_matrices) {
  const tensor::Tensor<int, 2, 3> x = {{
      {{7, 8, 9}},
      {{17, 18, 19}},
  }};

  ASSERT_THAT(as_vector(&x).rows(), Eq(6));
  ASSERT_THAT(as_vector(&x).cols(), Eq(1));
  EXPECT_THAT(as_vector(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_vector(&x)(1, 0), Eq(8));
  EXPECT_THAT(as_vector(&x)(2, 0), Eq(9));
  EXPECT_THAT(as_vector(&x)(3, 0), Eq(17));
  EXPECT_THAT(as_vector(&x)(4, 0), Eq(18));
  EXPECT_THAT(as_vector(&x)(5, 0), Eq(19));
}

TEST_F(as_vector_Test, works_for_nonconst_matrices) {
  tensor::Tensor<int, 2, 3> x = {{
      {{7, 8, 9}},
      {{17, 18, 19}},
  }};

  ASSERT_THAT(as_vector(&x).rows(), Eq(6));
  ASSERT_THAT(as_vector(&x).cols(), Eq(1));
  EXPECT_THAT(as_vector(&x)(0, 0), Eq(7));
  EXPECT_THAT(as_vector(&x)(1, 0), Eq(8));
  EXPECT_THAT(as_vector(&x)(2, 0), Eq(9));
  EXPECT_THAT(as_vector(&x)(3, 0), Eq(17));
  EXPECT_THAT(as_vector(&x)(4, 0), Eq(18));
  EXPECT_THAT(as_vector(&x)(5, 0), Eq(19));
}

TEST_F(as_vector_Test, modifying_nonconst_matrix_works) {
  tensor::Tensor<int, 2, 3> x = {{
      {{7, 8, 9}},
      {{17, 18, 19}},
  }};

  as_vector (&x)(0, 0) = 3;

  EXPECT_THAT(x.at(0).at(0), Eq(3));
}

} // namespace
} // namespace tensor
} // namespace elements
} // namespace ae108