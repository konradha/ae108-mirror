// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include <vector>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Returns a vector of size `N`, where `N` is the number of vertices in
 * the mesh. Element `i` of this vector is the index of the first degree of
 * freedom associated with vertex `i` in a distributed/global vector.
 */
template <class Policy>
std::vector<typename cpppetsc::Mesh<Policy>::size_type>
vertexDataOffsets(const cpppetsc::Mesh<Policy> &mesh);

extern template std::vector<
    typename cpppetsc::Mesh<SequentialComputePolicy>::size_type>
vertexDataOffsets(const cpppetsc::Mesh<SequentialComputePolicy> &);
extern template std::vector<
    typename cpppetsc::Mesh<ParallelComputePolicy>::size_type>
vertexDataOffsets(const cpppetsc::Mesh<ParallelComputePolicy> &);

} // namespace cpppetsc
} // namespace ae108

namespace ae108 {
namespace cpppetsc {

template <class Policy>
std::vector<typename cpppetsc::Mesh<Policy>::size_type>
vertexDataOffsets(const cpppetsc::Mesh<Policy> &mesh) {
  using Mesh = cpppetsc::Mesh<Policy>;
  using size_type = typename Mesh::size_type;

  auto offsets = std::vector<size_type>(mesh.totalNumberOfVertices());
  for (auto &&vertex : mesh.localVertices()) {
    offsets[vertex.index()] = vertex.globalDofLineRange().first;
  }
  Policy::handleError(MPI_Allreduce(MPI_IN_PLACE, offsets.data(),
                                    offsets.size(), MPIU_INT, MPIU_MAX,
                                    Policy::communicator()));
  return offsets;
}

} // namespace cpppetsc
} // namespace ae108