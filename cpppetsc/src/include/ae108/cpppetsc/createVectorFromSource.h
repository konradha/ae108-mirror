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

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include <functional>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Creates a global vector with `dofs` degrees of freedom per vertex in
 * the provided in the mesh. The vector is filled with data copied (not added)
 * from the source.
 *
 * @param mesh The mesh that defines the geometry, e.g. the number of vertices.
 *
 * @param dofs The desired number of degrees of freedom per vertex.
 *
 * @param writeVertexData For every local vertex, the function
 * is called with the vertex index and a buffer of size `dofs`. The function
 * shall write the data for the given vertex to the buffer.
 */
template <class Policy>
distributed<typename Mesh<Policy>::vector_type>
createVectorFromSource(const Mesh<Policy> &mesh,
                       const typename Mesh<Policy>::size_type dofs,
                       std::function<void(typename Mesh<Policy>::size_type,
                                          typename Mesh<Policy>::value_type *)>
                           writeVertexData);

extern template distributed<typename Mesh<SequentialComputePolicy>::vector_type>
createVectorFromSource(
    const Mesh<SequentialComputePolicy> &,
    Mesh<SequentialComputePolicy>::size_type,
    std::function<void(typename Mesh<SequentialComputePolicy>::size_type,
                       typename Mesh<SequentialComputePolicy>::value_type *)>);

extern template distributed<typename Mesh<ParallelComputePolicy>::vector_type>
createVectorFromSource(
    const Mesh<ParallelComputePolicy> &, Mesh<ParallelComputePolicy>::size_type,
    std::function<void(typename Mesh<ParallelComputePolicy>::size_type,
                       typename Mesh<ParallelComputePolicy>::value_type *)>);

} // namespace cpppetsc
} // namespace ae108

#include <vector>

namespace ae108 {
namespace cpppetsc {

template <class Policy>
distributed<typename Mesh<Policy>::vector_type>
createVectorFromSource(const Mesh<Policy> &mesh,
                       const typename Mesh<Policy>::size_type dofs,
                       std::function<void(typename Mesh<Policy>::size_type,
                                          typename Mesh<Policy>::value_type *)>
                           writeVertexData) {
  using mesh_type = Mesh<Policy>;
  using vector_type = typename mesh_type::vector_type;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;

  const auto elementDofs = size_type{0};
  auto layout = mesh.cloneWithDofs(dofs, elementDofs);
  auto local = vector_type::fromLocalMesh(layout);

  std::vector<value_type> data(dofs);
  for (const auto &vertex : layout.localVertices()) {
    writeVertexData(vertex.index(), data.data());
    vertex.setVertexData(data, &local);
  }

  auto global = vector_type::fromGlobalMesh(layout);
  layout.copyToGlobalVector(local, &global);

  return global;
}

} // namespace cpppetsc
} // namespace ae108
