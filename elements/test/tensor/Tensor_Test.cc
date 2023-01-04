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
#include <gmock/gmock.h>

using testing::Test;

namespace ae108 {
namespace elements {
namespace tensor {
namespace {

struct Tensor_Test : Test {};

TEST_F(Tensor_Test, zero_element_tensor_works) {
  static_assert(std::is_same<Tensor<double, 0>, std::array<double, 0>>::value,
                "A zero tensor is an array of length 0.");
}

TEST_F(Tensor_Test, one_element_tensor_works) {
  static_assert(std::is_same<Tensor<double, 1>, std::array<double, 1>>::value,
                "A zero tensor is an array of length 1.");
}

TEST_F(Tensor_Test, two_dimensional_tensor_works) {
  static_assert(std::is_same<Tensor<double, 1, 2>,
                             std::array<std::array<double, 2>, 1>>::value,
                "A two-dimensional tensor is an array of arrays.");
}

TEST_F(Tensor_Test, three_dimensional_tensor_works) {
  static_assert(
      std::is_same<Tensor<double, 1, 2, 3>,
                   std::array<std::array<std::array<double, 3>, 2>, 1>>::value,
      "A two-dimensional tensor is an array of arrays.");
}

} // namespace
} // namespace tensor
} // namespace elements
} // namespace ae108