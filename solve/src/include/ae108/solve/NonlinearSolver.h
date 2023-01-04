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

#include "ae108/assembly/AssemblerTypeTraits.h"
#include "ae108/cpppetsc/MeshBoundaryCondition.h"
#include "ae108/cpppetsc/NonlinearSolver.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include <cassert>
#include <functional>
#include <vector>

namespace ae108 {
namespace solve {

template <class Assembler> class NonlinearSolver {
public:
  using policy_type = typename assembly::PolicyTypeTrait<Assembler>::type;
  using mesh_type = typename assembly::MeshTypeTrait<Assembler>::type;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;
  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;

  /**
   * @param mesh A valid nonzero pointer to a mesh_type instance.
   */
  explicit NonlinearSolver(const mesh_type *mesh);

  using BoundaryConditionContainer =
      std::vector<cpppetsc::MeshBoundaryCondition<mesh_type>>;

  /**
   * @brief Calls computeSolution with the default assembler calls (ie. calls
   * assembleForceVector or assembleStiffnessMatrix). Passes on the rest of the
   * arguments.
   *
   * @param boundaryConditions The local essential boundary conditions to apply.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @para assembler Valid nonzero pointer.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(const BoundaryConditionContainer &boundaryConditions,
                  cpppetsc::distributed<vector_type> initialGuess,
                  const double time, const Assembler *const assembler) const;

  using LocalForceVectorAssembler =
      std::function<void(const cpppetsc::local<vector_type> &, double,
                         cpppetsc::local<vector_type> *)>;

  using LocalStiffnessMatrixAssembler = std::function<void(
      const cpppetsc::local<vector_type> &, double, matrix_type *)>;

  /**
   * @brief Calls computeSolution wrapping up the local assembler calls into
   * distributed assembler calls. Passes on the rest of the arguments.
   *
   * @param boundaryConditions The local essential boundary conditions to apply.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the local force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(const BoundaryConditionContainer &boundaryConditions,
                  cpppetsc::distributed<vector_type> initialGuess,
                  const double time,
                  LocalForceVectorAssembler assembleForceVector,
                  LocalStiffnessMatrixAssembler assemblerStiffnessMatrix) const;

  using DistributedForceVectorAssembler =
      std::function<void(const cpppetsc::distributed<vector_type> &, double,
                         cpppetsc::distributed<vector_type> *)>;

  using DistributedStiffnessMatrixAssembler = std::function<void(
      const cpppetsc::distributed<vector_type> &, double, matrix_type *)>;

  /**
   * @brief The goal of this function is to solve E(x) = min under the condition
   * that all the essential boundary conditions x_i = a_i are met.
   *
   * @remark Applying a nonlinear solver, we actually solve
   * F(x) := grad E(x) = 0 on the subspace where the boundary conditions are
   * satisfied. We turn this problem into a problem for all x by defining
   * f_i(x) := F_i(x) for subspace coordinates i and f_i(x) = x_i - a_i else.
   * The equation to solve therefore is:
   * f(x) = 0.
   *
   * @remark Does not update internal variables.
   *
   * @param boundaryConditions The local essential boundary conditions to apply.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the distributed force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type> computeSolution(
      const BoundaryConditionContainer &boundaryConditions,
      cpppetsc::distributed<vector_type> initialGuess, const double time,
      DistributedForceVectorAssembler assembleForceVector,
      DistributedStiffnessMatrixAssembler assemblerStiffnessMatrix) const;

private:
  const mesh_type *_mesh;
};
} // namespace solve
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace solve {
namespace detail {

template <class BoundaryConditionContainer>
std::vector<decltype(BoundaryConditionContainer::value_type::dof)>
extractGlobalIndices(const BoundaryConditionContainer &boundaryConditions) {
  using size_type = decltype(BoundaryConditionContainer::value_type::dof);
  std::vector<size_type> indices;
  indices.reserve(boundaryConditions.size());
  for (const auto &boundaryCondition : boundaryConditions) {
    const auto rowIndex = boundaryCondition.node.globalDofLineRange().first;
    const bool isLocallyOwned(rowIndex >= 0);
    if (isLocallyOwned) {
      indices.push_back(rowIndex + boundaryCondition.dof);
    }
  }
  return indices;
}

template <class BoundaryConditionContainer, class matrix_type>
void applyBoundaryConditionsToMatrix(
    const BoundaryConditionContainer &boundaryConditions,
    matrix_type *const matrix) {
  const auto indices = extractGlobalIndices(boundaryConditions);
  matrix->replaceRowsByEye(indices);
}

template <class vector_type, class BoundaryConditionContainer>
void applyBoundaryConditionsToForces(
    const BoundaryConditionContainer &boundaryConditions,
    const cpppetsc::distributed<vector_type> &input,
    cpppetsc::distributed<vector_type> *const output) {
  using size_type = decltype(BoundaryConditionContainer::value_type::dof);
  using value_type = decltype(BoundaryConditionContainer::value_type::value);

  auto indices = std::vector<size_type>();
  indices.reserve(boundaryConditions.size());
  auto values = std::vector<value_type>();
  values.reserve(boundaryConditions.size());

  for (const auto &boundaryCondition : boundaryConditions) {
    const auto globalRowIndex =
        boundaryCondition.node.globalDofLineRange().first;

    if (globalRowIndex >= 0) {
      indices.push_back(globalRowIndex + boundaryCondition.dof);
      values.push_back(input(globalRowIndex + boundaryCondition.dof) -
                       boundaryCondition.value);
    }
  }

  output->unwrap().replace().elements(indices, values);
}
} // namespace detail

template <class Assembler>
NonlinearSolver<Assembler>::NonlinearSolver(const mesh_type *mesh)
    : _mesh(mesh) {}

template <class Assembler>
cpppetsc::distributed<typename NonlinearSolver<Assembler>::vector_type>
NonlinearSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    const Assembler *const assembler) const {
  assert(assembler);
  return computeSolution(
      boundaryConditions, std::move(initialGuess), time,
      LocalForceVectorAssembler(
          [assembler](const cpppetsc::local<vector_type> &localDisplacements,
                      const double time,
                      cpppetsc::local<vector_type> *const localForces) {
            assembler->assembleForceVector(localDisplacements, time,
                                           localForces);
          }),
      LocalStiffnessMatrixAssembler(
          [assembler](const cpppetsc::local<vector_type> &localDisplacements,
                      const double time, matrix_type *const output) {
            assembler->assembleStiffnessMatrix(localDisplacements, time,
                                               output);
          }));
}

template <class Assembler>
cpppetsc::distributed<typename NonlinearSolver<Assembler>::vector_type>
NonlinearSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    LocalForceVectorAssembler assembleForceVector,
    LocalStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  const auto &mesh = *_mesh;

  auto localDisplacements = vector_type::fromLocalMesh(mesh);
  auto localForces = vector_type::fromLocalMesh(mesh);

  const auto distributedForceVectorAssembler =
      [&assembleForceVector, &mesh, &localDisplacements, &localForces](
          const cpppetsc::distributed<vector_type> &input, const double time,
          cpppetsc::distributed<vector_type> *const output) {
        mesh.copyToLocalVector(input, &localDisplacements);

        localForces.unwrap().setZero();
        assembleForceVector(localDisplacements, time, &localForces);

        mesh.addToGlobalVector(localForces, output);
      };

  const auto distributedStiffnessMatrixAssembler =
      [&assembleStiffnessMatrix, &mesh,
       &localDisplacements](const cpppetsc::distributed<vector_type> &input,
                            const double time, matrix_type *const output) {
        mesh.copyToLocalVector(input, &localDisplacements);

        assembleStiffnessMatrix(localDisplacements, time, output);
        output->finalize();
      };

  return computeSolution(boundaryConditions, std::move(initialGuess), time,
                         distributedForceVectorAssembler,
                         distributedStiffnessMatrixAssembler);
}

template <class Assembler>
cpppetsc::distributed<typename NonlinearSolver<Assembler>::vector_type>
NonlinearSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  auto function = [time, &assembleForceVector, &boundaryConditions](
                      const cpppetsc::distributed<vector_type> &input,
                      cpppetsc::distributed<vector_type> *const output) {
    output->unwrap().setZero();
    assembleForceVector(input, time, output);
    detail::applyBoundaryConditionsToForces<vector_type>(boundaryConditions,
                                                         input, output);
  };

  auto jacobian = [time, &assembleStiffnessMatrix, &boundaryConditions](
                      const cpppetsc::distributed<vector_type> &input,
                      matrix_type *const output) {
    output->setZero();
    assembleStiffnessMatrix(input, time, output);
    detail::applyBoundaryConditionsToMatrix(boundaryConditions, output);
  };

  const cpppetsc::NonlinearSolver<policy_type> solver{
      matrix_type::fromMesh(*_mesh), vector_type::fromGlobalMesh(*_mesh)};
  return solver.solve(std::move(function), std::move(jacobian),
                      std::move(initialGuess));
}

} // namespace solve
} // namespace ae108
