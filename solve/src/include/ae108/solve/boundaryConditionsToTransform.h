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

#include "ae108/cpppetsc/GeneralizedMeshBoundaryCondition.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/createRhsTransform.h"
#include "ae108/cpppetsc/createTransformOutput.h"
#include "ae108/cpppetsc/vertexDataOffsets.h"
#include "ae108/solve/AffineTransform.h"
#include "ae108/solve/InconsistentBoundaryConditionsException.h"
#include "ae108/solve/InvalidVertexException.h"
#include <algorithm>
#include <utility>
#include <vector>

namespace ae108 {
namespace solve {

/**
 * @brief Converts a vector of `GeneralizedMeshBoundaryCondition`s to an affine
 * transform that encodes the boundary conditions. More precisely, the rank of
 * the transform satisfies the boundary conditions.
 *
 * @throw InconsistentBoundaryConditionsException if a degree of degree of
 * freedom appears as both a source and a target.
 *
 * @throw InvalidVertexException if the boundary conditions
 * contain a vertex that does not exist.
 */
template <class Policy>
AffineTransform<Policy> boundaryConditionsToTransform(
    const std::vector<
        cpppetsc::GeneralizedMeshBoundaryCondition<cpppetsc::Mesh<Policy>>>
        &boundaryConditions,
    const cpppetsc::Mesh<Policy> &mesh);

extern template AffineTransform<cpppetsc::SequentialComputePolicy>
boundaryConditionsToTransform(
    const std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<
        cpppetsc::Mesh<cpppetsc::SequentialComputePolicy>>> &,
    const cpppetsc::Mesh<cpppetsc::SequentialComputePolicy> &);

extern template AffineTransform<cpppetsc::ParallelComputePolicy>
boundaryConditionsToTransform(
    const std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<
        cpppetsc::Mesh<cpppetsc::ParallelComputePolicy>>> &,
    const cpppetsc::Mesh<cpppetsc::ParallelComputePolicy> &);

} // namespace solve
} // namespace ae108

namespace ae108 {
namespace solve {

template <class Policy>
AffineTransform<Policy> boundaryConditionsToTransform(
    const std::vector<
        cpppetsc::GeneralizedMeshBoundaryCondition<cpppetsc::Mesh<Policy>>>
        &boundaryConditions,
    const cpppetsc::Mesh<Policy> &mesh) {
  using Mesh = cpppetsc::Mesh<Policy>;
  using BoundaryCondition = cpppetsc::GeneralizedMeshBoundaryCondition<Mesh>;
  using size_type = typename Mesh::size_type;
  using value_type = typename Mesh::value_type;

  const auto vertexDataOffsets = cpppetsc::vertexDataOffsets(mesh);
  const auto itemToRow = [&](const typename BoundaryCondition::Item &item) {
    try {
      return vertexDataOffsets.at(item.vertex) + item.dof;
    } catch (std::out_of_range &) {
      throw InvalidVertexException{};
    }
  };

  const auto matrix = Mesh::matrix_type::fromMesh(mesh);

  const auto eliminatedRows = [&](const size_type numberOfRows) {
    std::vector<size_type> rows(numberOfRows);
    for (auto &&condition : boundaryConditions) {
      rows[itemToRow(condition.target)] = size_type{1};
    }
    Policy::handleError(MPI_Allreduce(MPI_IN_PLACE, rows.data(), rows.size(),
                                      MPIU_INT, MPIU_MAX,
                                      Policy::communicator()));
    return rows;
  }(matrix.size().first);

  const auto offsets = [&]() {
    auto offsets = std::vector<size_type>(eliminatedRows.size());
    std::partial_sum(eliminatedRows.begin(), eliminatedRows.end(),
                     offsets.begin());
    return offsets;
  }();

  auto transform = cpppetsc::createRhsTransform(matrix, matrix.size().second -
                                                            offsets.back());
  auto shift = cpppetsc::createTransformOutput(transform);

  {
    const auto localRange = matrix.localRowRange();
    auto matrixReplacer = transform.assemblyView().replace();
    for (auto row = localRange.first; row < localRange.second; ++row) {
      if (!eliminatedRows[row]) {
        const auto col = row;
        matrixReplacer(row, col - offsets[col]) = value_type{1.};
      }
    }

    auto vectorReplacer = shift.unwrap().replace();
    for (auto &&condition : boundaryConditions) {
      const auto row = itemToRow(condition.target);
      for (auto &&source : condition.source) {
        const auto col = itemToRow(source.item);
        if (eliminatedRows[col]) {
          throw InconsistentBoundaryConditionsException{};
        }
        matrixReplacer(row, col - offsets[col]) = source.factor;
      }
      vectorReplacer(row) = condition.offset;
    }
  }

  return AffineTransform<Policy>{std::move(transform), std::move(shift)};
}
} // namespace solve
} // namespace ae108