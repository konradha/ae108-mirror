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

TEST_F(Tensor_Test, tensor_of_float_works) {
  static_assert(std::is_same<Tensor<float, 1>, std::array<float, 1>>::value,
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