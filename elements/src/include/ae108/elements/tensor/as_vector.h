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