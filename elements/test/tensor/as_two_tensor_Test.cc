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
#include "ae108/elements/tensor/as_two_tensor.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace elements {
namespace tensor {
namespace {

struct as_two_tensor_Test : Test {
  using Data = std::array<std::array<std::array<std::array<int, 4>, 3>, 2>, 1>;
  Data data = {{
      {{
          {{
              {{1, 2, 3, 4}},
              {{5, 6, 7, 8}},
              {{9, 10, 11, 12}},
          }},
          {{
              {{11, 12, 13, 14}},
              {{15, 16, 17, 18}},
              {{19, 110, 111, 112}},
          }},
      }},
  }};
  const Data const_data = data;
};

TEST_F(as_two_tensor_Test, has_correct_dimensions) {
  EXPECT_THAT(as_two_tensor(&data).rows(), Eq(2 * 1));
  EXPECT_THAT(as_two_tensor(&data).cols(), Eq(4 * 3));
}

TEST_F(as_two_tensor_Test, const_has_correct_dimensions) {
  EXPECT_THAT(as_two_tensor(&const_data).rows(),
              Eq(as_two_tensor(&data).rows()));
  EXPECT_THAT(as_two_tensor(&const_data).cols(),
              Eq(as_two_tensor(&data).cols()));
}

TEST_F(as_two_tensor_Test, has_correct_entries) {
  const auto matrix = as_two_tensor(&data);
  EXPECT_THAT(matrix(0 * 1 + 0, 0 * 4 + 0), Eq(data[0][0][0][0]));
  EXPECT_THAT(matrix(0 * 1 + 0, 0 * 4 + 1), Eq(data[0][0][0][1]));
  EXPECT_THAT(matrix(0 * 1 + 0, 0 * 4 + 2), Eq(data[0][0][0][2]));
  EXPECT_THAT(matrix(0 * 1 + 0, 0 * 4 + 3), Eq(data[0][0][0][3]));
  EXPECT_THAT(matrix(0 * 1 + 0, 1 * 4 + 0), Eq(data[0][0][1][0]));
  EXPECT_THAT(matrix(0 * 1 + 0, 1 * 4 + 1), Eq(data[0][0][1][1]));
  EXPECT_THAT(matrix(0 * 1 + 0, 1 * 4 + 2), Eq(data[0][0][1][2]));
  EXPECT_THAT(matrix(0 * 1 + 0, 1 * 4 + 3), Eq(data[0][0][1][3]));
  EXPECT_THAT(matrix(0 * 1 + 0, 2 * 4 + 0), Eq(data[0][0][2][0]));
  EXPECT_THAT(matrix(0 * 1 + 0, 2 * 4 + 1), Eq(data[0][0][2][1]));
  EXPECT_THAT(matrix(0 * 1 + 0, 2 * 4 + 2), Eq(data[0][0][2][2]));
  EXPECT_THAT(matrix(0 * 1 + 0, 2 * 4 + 3), Eq(data[0][0][2][3]));
  EXPECT_THAT(matrix(0 * 1 + 1, 0 * 4 + 0), Eq(data[0][1][0][0]));
  EXPECT_THAT(matrix(0 * 1 + 1, 0 * 4 + 1), Eq(data[0][1][0][1]));
  EXPECT_THAT(matrix(0 * 1 + 1, 0 * 4 + 2), Eq(data[0][1][0][2]));
  EXPECT_THAT(matrix(0 * 1 + 1, 0 * 4 + 3), Eq(data[0][1][0][3]));
  EXPECT_THAT(matrix(0 * 1 + 1, 1 * 4 + 0), Eq(data[0][1][1][0]));
  EXPECT_THAT(matrix(0 * 1 + 1, 1 * 4 + 1), Eq(data[0][1][1][1]));
  EXPECT_THAT(matrix(0 * 1 + 1, 1 * 4 + 2), Eq(data[0][1][1][2]));
  EXPECT_THAT(matrix(0 * 1 + 1, 1 * 4 + 3), Eq(data[0][1][1][3]));
  EXPECT_THAT(matrix(0 * 1 + 1, 2 * 4 + 0), Eq(data[0][1][2][0]));
  EXPECT_THAT(matrix(0 * 1 + 1, 2 * 4 + 1), Eq(data[0][1][2][1]));
  EXPECT_THAT(matrix(0 * 1 + 1, 2 * 4 + 2), Eq(data[0][1][2][2]));
  EXPECT_THAT(matrix(0 * 1 + 1, 2 * 4 + 3), Eq(data[0][1][2][3]));
}

TEST_F(as_two_tensor_Test, const_has_correct_entries) {
  EXPECT_THAT(as_two_tensor(&const_data), Eq(as_two_tensor(&data)));
}

TEST_F(as_two_tensor_Test, data_can_be_written) {
  const auto value = 777;

  auto matrix = as_two_tensor(&data);
  matrix(0, 0) = value;

  EXPECT_THAT(data[0][0][0][0], Eq(value));
}

} // namespace
} // namespace tensor
} // namespace elements
} // namespace ae108