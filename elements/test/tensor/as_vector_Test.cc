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

} // namespace
} // namespace tensor
} // namespace elements
} // namespace ae108