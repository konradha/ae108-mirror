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
