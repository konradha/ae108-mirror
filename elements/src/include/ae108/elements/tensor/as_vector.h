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
 * @brief Interprets the input as a vector and returns an Eigen
 * vector.
 */
template <class ValueType_, std::size_t Rows_>
Eigen::Map<Eigen::Matrix<ValueType_, Rows_, 1>>
as_vector(std::array<ValueType_, Rows_> *data) {
  assert(data);
  return Eigen::Map<Eigen::Matrix<ValueType_, Rows_, 1>>(data->data());
}

/**
 * @brief Interprets the input as a vector and returns a constant Eigen
 * vector.
 */
template <class ValueType_, std::size_t Rows_>
Eigen::Map<const Eigen::Matrix<ValueType_, Rows_, 1>>
as_vector(const std::array<ValueType_, Rows_> *data) {
  assert(data);
  return Eigen::Map<const Eigen::Matrix<ValueType_, Rows_, 1>>(data->data());
}

/**
 * @brief Interprets the input as a vector of columns stacked on top of each
 * other and returns an Eigen vector.
 */
template <class ValueType_, std::size_t Rows_, std::size_t Cols_>
Eigen::Map<Eigen::Matrix<ValueType_, Rows_ * Cols_, 1>>
as_vector(std::array<std::array<ValueType_, Cols_>, Rows_> *data) {
  static_assert(Cols_ * Rows_ == 0 || sizeof(decltype(*data)) ==
                                          Cols_ * Rows_ * sizeof(ValueType_),
                "Only contiguous data is supported.");
  assert(data);
  return Eigen::Map<Eigen::Matrix<ValueType_, Rows_ * Cols_, 1>>(
      Rows_ > 0 ? data->front().data() : nullptr);
}

/**
 * @brief Interprets the input as a vector of columns stacked on top of each
 * other and returns a constant Eigen vector.
 */
template <class ValueType_, std::size_t Rows_, std::size_t Cols_>
Eigen::Map<const Eigen::Matrix<ValueType_, Rows_ * Cols_, 1>>
as_vector(const std::array<std::array<ValueType_, Cols_>, Rows_> *data) {
  static_assert(Cols_ * Rows_ == 0 || sizeof(decltype(*data)) ==
                                          Cols_ * Rows_ * sizeof(ValueType_),
                "Only contiguous data is supported.");
  assert(data);
  return Eigen::Map<const Eigen::Matrix<ValueType_, Rows_ * Cols_, 1>>(
      Rows_ > 0 ? data->front().data() : nullptr);
}

} // namespace tensor
} // namespace elements
} // namespace ae108