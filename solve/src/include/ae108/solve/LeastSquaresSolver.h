// © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
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

#ifndef AE108_PETSC_COMPLEX

#include "ae108/assembly/AssemblerTypeTraits.h"
#include "ae108/cpppetsc/LeastSquaresSolver.h"
#include "ae108/cpppetsc/MeshBoundaryCondition.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/clone.h"
#include "ae108/solve/AffineTransform.h"
#include <cassert>
#include <functional>
#include <vector>

namespace ae108 {
namespace solve {

template <class Assembler> class LeastSquaresSolver {
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
  explicit LeastSquaresSolver(const mesh_type *mesh);

  using BoundaryConditionContainer = AffineTransform<policy_type>;

  /**
   * @brief Calls computeSolution with the default assembler calls (ie. calls
   * assembleForceVector or assembleStiffnessMatrix). Passes on the rest of the
   * arguments.
   *
   * @param transform The system of linear equations.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @para assembler Valid nonzero pointer.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(const BoundaryConditionContainer &transform,
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
   * @param transform The system of linear equations.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the local force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(const BoundaryConditionContainer &transform,
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
   * @brief The goal of this function is to find the minimizer x of
   * E(x) -> min, where the domain is defined by boundary conditions
   * specified by linear equations.
   *
   * More precisely, this function solves the equation E' = F = 0. Denote the
   * system of linear equations that corresponds to the boundary conditions by
   * Ax + b = 0. Then this solver minimizes |F(x)|^2 + |Ax + b|^2.
   *
   * @param transform The system of linear equations.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the distributed force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type> computeSolution(
      const BoundaryConditionContainer &transform,
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

#include "ae108/cpppetsc/createTransformInput.h"

namespace ae108 {
namespace solve {

template <class Assembler>
LeastSquaresSolver<Assembler>::LeastSquaresSolver(const mesh_type *mesh)
    : _mesh(mesh) {}

template <class Assembler>
cpppetsc::distributed<typename LeastSquaresSolver<Assembler>::vector_type>
LeastSquaresSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &transform,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    const Assembler *const assembler) const {
  assert(assembler);
  return computeSolution(
      transform, std::move(initialGuess), time,
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
cpppetsc::distributed<typename LeastSquaresSolver<Assembler>::vector_type>
LeastSquaresSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &transform,
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

  return computeSolution(transform, std::move(initialGuess), time,
                         distributedForceVectorAssembler,
                         distributedStiffnessMatrixAssembler);
}

template <class Assembler>
cpppetsc::distributed<typename LeastSquaresSolver<Assembler>::vector_type>
LeastSquaresSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &transform,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  auto residual = clone(transform.vector);

  auto fullForces = vector_type::fromGlobalMesh(*_mesh);
  auto nestedVector = [&]() {
    auto nestedVector = Vec{};
    auto vectors = std::array<Vec, 2>{
        {fullForces.unwrap().data(), residual.unwrap().data()}};
    policy_type::handleError(VecCreateNest(policy_type::communicator(), 2,
                                           nullptr, vectors.data(),
                                           &nestedVector));
    return cpppetsc::distributed<vector_type>(
        cpppetsc::makeUniqueEntity<policy_type>(nestedVector));
  }();

  auto function = [&](const cpppetsc::distributed<vector_type> &input,
                      cpppetsc::distributed<vector_type> *const) {
    fullForces.unwrap().setZero();
    assembleForceVector(input, time, &fullForces);
    apply(transform, input, &residual);
  };

  auto fullStiffnessMatrix = matrix_type::fromMesh(*_mesh);
  auto nestedMatrix = [&]() {
    auto nestedMatrix = Mat{};
    const auto matrices = std::array<Mat, 2>{
        {fullStiffnessMatrix.data(), transform.matrix.data()}};
    policy_type::handleError(MatCreateNest(policy_type::communicator(), 2,
                                           nullptr, 1, nullptr, matrices.data(),
                                           &nestedMatrix));
    return matrix_type(cpppetsc::makeUniqueEntity<policy_type>(nestedMatrix));
  }();

  auto jacobian = [&](const cpppetsc::distributed<vector_type> &input,
                      matrix_type *const) {
    fullStiffnessMatrix.setZero();
    assembleStiffnessMatrix(input, time, &fullStiffnessMatrix);
  };

  const cpppetsc::LeastSquaresSolver<policy_type> solver{
      std::move(nestedMatrix), std::move(nestedVector)};

  return solver.solve(std::move(function), std::move(jacobian),
                      std::move(initialGuess));
}

} // namespace solve
} // namespace ae108

#endif