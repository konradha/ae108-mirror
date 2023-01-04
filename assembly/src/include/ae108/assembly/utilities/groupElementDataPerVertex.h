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

#include "ae108/assembly/utilities/deserialize.h"
#include "ae108/assembly/utilities/resizeIfPossible.h"
#include "ae108/assembly/utilities/serialize.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include <cassert>
#include <numeric>
#include <vector>

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief Copies the data from the input vector to the output, grouped by
 * vertex.
 *
 * @param to Valid nonzero pointer.
 */
template <class PerVertexData, class MeshView, class VectorType>
void groupElementDataPerVertex(const MeshView &meshView,
                               const cpppetsc::local<VectorType> &from,
                               PerVertexData *const to) {
  assert(to);

  std::vector<typename VectorType::value_type> buffer;
  meshView.copyElementData(from, &buffer);

  resizeIfPossible(to, meshView.numberOfVertices());
  deserializeRange(buffer.begin(), buffer.end(), to->begin());
}

/**
 * @brief Copies the data from the input vector to the output, ungrouped by
 * vertex.
 *
 * @param to Valid nonzero pointer.
 */
template <class PerVertexData, class MeshView, class VectorType>
void ungroupElementDataPerVertex(const MeshView &meshView,
                                 const PerVertexData &from,
                                 cpppetsc::local<VectorType> *const to) {
  assert(to);

  using BufferType = std::vector<typename VectorType::value_type>;
  BufferType buffer(std::accumulate(
      from.begin(), from.end(), typename BufferType::size_type{},
      [](const typename BufferType::size_type sum,
         const typename PerVertexData::value_type &value) {
        return sum + value.size();
      }));

  serializeRange(from.begin(), from.end(), buffer.begin());
  meshView.addElementData(buffer, to);
}
} // namespace utilities
} // namespace assembly
} // namespace ae108
