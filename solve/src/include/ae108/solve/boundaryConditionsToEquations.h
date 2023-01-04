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

#include "ae108/cpppetsc/GeneralizedMeshBoundaryCondition.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/createLhsTransform.h"
#include "ae108/cpppetsc/createTransformOutput.h"
#include "ae108/cpppetsc/vertexDataOffsets.h"
#include "ae108/solve/AffineTransform.h"
#include "ae108/solve/InvalidVertexException.h"
#include <algorithm>
#include <utility>
#include <vector>

namespace ae108 {
namespace solve {

/**
 * @brief Converts a vector of `GeneralizedMeshBoundaryCondition`s to a matrix
 * `A` and a vector `b` that encodes the boundary conditions. The boundary
 * conditions are satisfied for `x` if and only if `Ax + b = 0`.
 *
 * @throw InvalidVertexException if the boundary conditions
 * contain a vertex that does not exist.
 */
template <class Policy>
AffineTransform<Policy> boundaryConditionsToEquations(
    const std::vector<
        cpppetsc::GeneralizedMeshBoundaryCondition<cpppetsc::Mesh<Policy>>>
        &boundaryConditions,
    const cpppetsc::Mesh<Policy> &mesh);

extern template AffineTransform<cpppetsc::SequentialComputePolicy>
boundaryConditionsToEquations(
    const std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<
        cpppetsc::Mesh<cpppetsc::SequentialComputePolicy>>> &,
    const cpppetsc::Mesh<cpppetsc::SequentialComputePolicy> &);

extern template AffineTransform<cpppetsc::ParallelComputePolicy>
boundaryConditionsToEquations(
    const std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<
        cpppetsc::Mesh<cpppetsc::ParallelComputePolicy>>> &,
    const cpppetsc::Mesh<cpppetsc::ParallelComputePolicy> &);

} // namespace solve
} // namespace ae108

namespace ae108 {
namespace solve {

template <class Policy>
AffineTransform<Policy> boundaryConditionsToEquations(
    const std::vector<
        cpppetsc::GeneralizedMeshBoundaryCondition<cpppetsc::Mesh<Policy>>>
        &boundaryConditions,
    const cpppetsc::Mesh<Policy> &mesh) {
  using Mesh = cpppetsc::Mesh<Policy>;
  using BoundaryCondition = cpppetsc::GeneralizedMeshBoundaryCondition<Mesh>;
  using size_type = typename Mesh::size_type;
  using value_type = typename Mesh::value_type;
  using matrix_type = typename Mesh::matrix_type;

  const auto vertexDataOffsets = cpppetsc::vertexDataOffsets(mesh);
  const auto itemToCol = [&](const typename BoundaryCondition::Item &item) {
    try {
      return vertexDataOffsets.at(item.vertex) + item.dof;
    } catch (std::out_of_range &) {
      throw InvalidVertexException{};
    }
  };

  const auto matrix = matrix_type::fromMesh(mesh);
  auto transform = createLhsTransform(
      matrix, typename matrix_type::LocalRows{boundaryConditions.size()},
      typename matrix_type::GlobalRows{PETSC_DETERMINE});
  auto shift = cpppetsc::createTransformOutput(transform);

  {
    const auto nonZeroesPerRow =
        size_type{1} +
        (boundaryConditions.empty()
             ? size_type{0}
             : static_cast<size_type>(
                   std::max_element(boundaryConditions.begin(),
                                    boundaryConditions.end(),
                                    [](const BoundaryCondition &x,
                                       const BoundaryCondition &y) {
                                      return x.source.size() < y.source.size();
                                    })
                       ->source.size()));

    auto matrixReplacer =
        transform.preallocatedAssemblyView(nonZeroesPerRow).replace();
    auto vectorReplacer = shift.unwrap().replace();

    auto range = transform.localRowRange();
    for (auto &&condition : boundaryConditions) {
      assert(range.first < range.second);

      matrixReplacer(range.first, itemToCol(condition.target)) = value_type{1.};
      vectorReplacer(range.first) = value_type{-1.} * condition.offset;

      for (auto &&source : condition.source) {
        matrixReplacer(range.first, itemToCol(source.item)) =
            value_type{-1.} * source.factor;
      }

      range.first++;
    }
  }

  return AffineTransform<Policy>{std::move(transform), std::move(shift)};
}
} // namespace solve
} // namespace ae108