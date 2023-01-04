// © 2021 ETH Zurich, Mechanics and Materials Lab
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