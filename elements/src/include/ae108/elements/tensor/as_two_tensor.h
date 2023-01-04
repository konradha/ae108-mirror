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