// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief Deserialize a buffer to an array by copying the buffer to the array.
 *
 * @tparam InputIterator Forward iterator.
 * @param bufferBegin Data source iterator.
 * @param array Data target.
 *
 * @return The end of the used input buffer.
 */
template <class InputIterator, class ValueType, std::size_t ArraySize>
InputIterator deserialize(InputIterator bufferBegin,
                          std::array<ValueType, ArraySize> *const array) {
  assert(array);
  auto bufferEnd = bufferBegin;
  std::advance(bufferEnd, ArraySize);
  std::copy(bufferBegin, bufferEnd, array->begin());
  return bufferEnd;
}

/**
 * @brief Deserialize a buffer to a vector by copying the buffer to the vector.
 *
 * @tparam InputIterator Forward iterator.
 * @param bufferBegin Data source iterator.
 * @param vector Data target.
 *
 * @return The end of the used input buffer.
 */
template <class InputIterator, class ValueType>
InputIterator deserialize(InputIterator bufferBegin,
                          std::vector<ValueType> *const vector) {
  assert(vector);
  auto bufferEnd = bufferBegin;
  std::advance(bufferEnd, vector->size());
  std::copy(bufferBegin, bufferEnd, vector->begin());
  return bufferEnd;
}

/**
 * @brief Deserialize a buffer to an Eigen::Matrix by copying the buffer to the
 * matrix (row major format).
 *
 * @tparam InputIterator Forward iterator.
 * @param bufferBegin Data source iterator.
 * @param vector Data target.
 *
 * @return The end of the used input buffer.
 */
template <class InputIterator, class ValueType, int Rows, int Cols>
InputIterator deserialize(InputIterator bufferBegin,
                          Eigen::Matrix<ValueType, Rows, Cols> *vector) {
  assert(vector);
  std::array<ValueType, Rows * Cols> buffer;
  const auto result = deserialize(bufferBegin, &buffer);
  *vector = Eigen::Map<const Eigen::Matrix<
      ValueType, Rows, Cols, Cols == 1 ? Eigen::ColMajor : Eigen::RowMajor>>(
      buffer.data(), vector->rows(), vector->cols());
  return result;
}

/**
 * @brief Iterates through the input range and deserializes chunks of buffer.
 *
 * @tparam InputIterator Forward iterator.
 * @tparam OutputIterator Output iterator.
 * @param bufferBegin Data source begin iterator.
 * @param bufferEnd Data source end iterator.
 * @param outputBegin Data target begin iterator.
 *
 * @return The end of the used input buffer.
 */
template <class InputIterator, class OutputIterator>
InputIterator deserializeRange(InputIterator bufferBegin,
                               InputIterator bufferEnd,
                               OutputIterator outputBegin) {
  while (bufferBegin != bufferEnd) {
    bufferBegin = deserialize(bufferBegin, &*outputBegin);
    outputBegin++;
  }
  return bufferBegin;
}
} // namespace utilities
} // namespace assembly
} // namespace ae108
