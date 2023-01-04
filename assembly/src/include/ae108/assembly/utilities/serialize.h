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
  for (Eigen::Index i = 0; i < vector.rows(); ++i) {
    for (Eigen::Index j = 0; j < vector.cols(); ++j) {
      *bufferBegin = vector(i, j);
      bufferBegin++;
    }
  }
  return bufferBegin;
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
