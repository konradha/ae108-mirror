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

#pragma once

#include <Eigen/Core>
#include <cassert>
#include <type_traits>

namespace ae108 {
namespace elements {
namespace tensor {

/**
 * @brief Interprets the input 4-tensor as a row-major 2-tensor matrix.
 */
template <class ValueType_, std::size_t A_, std::size_t B_, std::size_t C_,
          std::size_t D_>
Eigen::Map<Eigen::Matrix<ValueType_, A_ * B_, C_ * D_,
                         (C_ * D_ == 1 && A_ * B_ != 1) ? Eigen::ColMajor
                                                        : Eigen::RowMajor>>
as_two_tensor(
    std::array<std::array<std::array<std::array<ValueType_, D_>, C_>, B_>, A_>
        *data) {
  static_assert(A_ * B_ * C_ * D_ == 0 ||
                    sizeof(decltype(*data)) ==
                        A_ * B_ * C_ * D_ * sizeof(ValueType_),
                "Only contiguous data is supported.");
  assert(data);
  return Eigen::Map<Eigen::Matrix<
      ValueType_, A_ * B_, C_ * D_,
      (C_ * D_ == 1 && A_ * B_ != 1) ? Eigen::ColMajor : Eigen::RowMajor>>(
      A_ > 0 && B_ > 0 && C_ > 0 ? data->front().front().front().data()
                                 : nullptr);
}

/**
 * @brief Interprets the input 4-tensor as a row-major 2-tensor matrix.
 */
template <class ValueType_, std::size_t A_, std::size_t B_, std::size_t C_,
          std::size_t D_>
Eigen::Map<const Eigen::Matrix<
    ValueType_, A_ * B_, C_ * D_,
    (C_ * D_ == 1 && A_ * B_ != 1) ? Eigen::ColMajor : Eigen::RowMajor>>
as_two_tensor(
    const std::array<std::array<std::array<std::array<ValueType_, D_>, C_>, B_>,
                     A_> *data) {
  static_assert(A_ * B_ * C_ * D_ == 0 ||
                    sizeof(decltype(*data)) ==
                        A_ * B_ * C_ * D_ * sizeof(ValueType_),
                "Only contiguous data is supported.");
  assert(data);
  return Eigen::Map<const Eigen::Matrix<
      ValueType_, A_ * B_, C_ * D_,
      (C_ * D_ == 1 && A_ * B_ != 1) ? Eigen::ColMajor : Eigen::RowMajor>>(
      A_ > 0 && B_ > 0 && C_ > 0 ? data->front().front().front().data()
                                 : nullptr);
}

} // namespace tensor
} // namespace elements
} // namespace ae108