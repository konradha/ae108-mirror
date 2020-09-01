// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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
#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief Serializes an array to a buffer by copying the array to the buffer.
 *
 * @tparam OutputIterator Forward iterator.
 * @param array Data source.
 * @param bufferBegin Data target begin iterator.
 *
 * @return The end of the output buffer.
 */
template <class OutputIterator, class ValueType, std::size_t ArraySize>
OutputIterator serialize(const std::array<ValueType, ArraySize> &array,
                         OutputIterator bufferBegin) {
  return std::copy(array.begin(), array.end(), bufferBegin);
}

/**
 * @brief Serializes a vector to a buffer by copying the vector to the buffer.
 *
 * @tparam OutputIterator Forward iterator.
 * @param vector Data source.
 * @param bufferBegin Data target begin iterator.
 *
 * @return The end of the output buffer.
 */
template <class OutputIterator, class ValueType>
OutputIterator serialize(const std::vector<ValueType> &vector,
                         OutputIterator bufferBegin) {
  return std::copy(vector.begin(), vector.end(), bufferBegin);
}

/**
 * @brief Serializes an Eigen::Matrix to a buffer by copying the vector to the
 * buffer (row major format).
 *
 * @tparam OutputIterator Forward iterator.
 * @param vector Data source.
 * @param bufferBegin Data target begin iterator.
 *
 * @return The end of the output buffer.
 */
template <class OutputIterator, class ValueType, int Rows, int Cols>
OutputIterator serialize(const Eigen::Matrix<ValueType, Rows, Cols> &vector,
                         OutputIterator bufferBegin) {
  std::array<ValueType, Rows * Cols> buffer;
  Eigen::Map<Eigen::Matrix<ValueType, Rows, Cols,
                           Cols == 1 ? Eigen::ColMajor : Eigen::RowMajor>>
      wrappedBuffer(buffer.data(), vector.rows(), vector.cols());
  wrappedBuffer = vector;
  return std::copy(buffer.begin(), buffer.end(), bufferBegin);
}

/**
 * @brief Iterates through the input range serializing every vector to the
 * output buffer.
 *
 * @tparam InputIterator Forward iterator.
 * @tparam OutputIterator Random access iterator.
 * @param vectorsBegin Data source begin iterator.
 * @param vectorsEnd Data source end iterator.
 * @param bufferBegin Data target begin iterator.
 *
 * @return The end of the used output buffer.
 */
template <class InputIterator, class OutputIterator>
OutputIterator serializeRange(InputIterator vectorsBegin,
                              InputIterator vectorsEnd,
                              OutputIterator bufferBegin) {
  while (vectorsBegin != vectorsEnd) {
    bufferBegin = serialize(*vectorsBegin, bufferBegin);
    ++vectorsBegin;
  }
  return bufferBegin;
}
} // namespace utilities
} // namespace assembly
} // namespace ae108
