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

#include "ae108/assembly/AssemblerTypeTraits.h"
#include "ae108/cpppetsc/GeneralizedMeshBoundaryCondition.h"
#include "ae108/cpppetsc/NonlinearSolver.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/asTransposedMatrix.h"
#include "ae108/solve/AffineTransform.h"
#include <functional>
#include <vector>

namespace ae108 {
namespace solve {

template <class Assembler> class GeneralizedNonlinearSolver {
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
  explicit GeneralizedNonlinearSolver(const mesh_type *mesh);

  using BoundaryConditionContainer =
      std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<mesh_type>>;

  /**
   * @brief Calls computeSolution with the default assembler calls (ie. calls
   * assembleForceVector or assembleStiffnessMatrix). Passes on the rest of the
   * arguments.
   *
   * @param boundaryConditions The generalized boundary conditions.
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
   * @param boundaryConditions The generalized boundary conditions.
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
   * @brief The goal of this function is to find the minimizer x of
   * E(x) -> min, where x satisfies boundary conditions specified by linear
   * equations.
   *
   * More precisely, this function uses an "infeasible start Newton
   * method"-based approach;  see S. Boyd, L. Vandenberghe, "Convex
   * Optimization", Cambridge University Press, 2019, pp. 531-534. PETSc's
   * nonlinear solver is used to solve the adapted system.
   *
   * This method is also widely known as Lagrange Multiplier Method;
   * see Cook, R. D., Malkus, D. S., Plesha, M. E. & Witt,
   * R. J. Concepts and Applications of Finite Element Analysis,
   * 4th Edition. (Wiley, 2001).
   *
   * @param boundaryConditions The generalized boundary conditions.
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

#include "ae108/cpppetsc/createTransformInput.h"
#include "ae108/cpppetsc/nestMatrices.h"
#include "ae108/cpppetsc/nestVectors.h"
#include "ae108/solve/boundaryConditionsToEquations.h"
#include <algorithm>
#include <cassert>

namespace ae108 {
namespace solve {

template <class Assembler>
GeneralizedNonlinearSolver<Assembler>::GeneralizedNonlinearSolver(
    const mesh_type *mesh)
    : _mesh(mesh) {}

template <class Assembler>
cpppetsc::distributed<
    typename GeneralizedNonlinearSolver<Assembler>::vector_type>
GeneralizedNonlinearSolver<Assembler>::computeSolution(
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
cpppetsc::distributed<
    typename GeneralizedNonlinearSolver<Assembler>::vector_type>
GeneralizedNonlinearSolver<Assembler>::computeSolution(
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
cpppetsc::distributed<
    typename GeneralizedNonlinearSolver<Assembler>::vector_type>
GeneralizedNonlinearSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  auto transform =
      solve::boundaryConditionsToEquations(boundaryConditions, *_mesh);
  auto transposed = asTransposedMatrix(&transform.matrix);

  auto fullForces = vector_type::fromGlobalMesh(*_mesh);
  auto dual = createTransformInput(transposed);
  auto function = [&](const cpppetsc::distributed<vector_type> &input,
                      cpppetsc::distributed<vector_type> *const) {
    fullForces.unwrap().setZero();
    assembleForceVector(input, time, &fullForces);

    apply(transform, initialGuess, &dual);
  };

  auto fullStiffnessMatrix = matrix_type::fromMesh(*_mesh);
  auto jacobian = [&](const cpppetsc::distributed<vector_type> &input,
                      matrix_type *const) {
    fullStiffnessMatrix.setZero();
    assembleStiffnessMatrix(input, time, &fullStiffnessMatrix);
  };

  auto zeroes = matrix_type(
      typename matrix_type::LocalRows{transform.matrix.localSize().first},
      typename matrix_type::LocalCols{transform.matrix.localSize().first},
      typename matrix_type::GlobalRows{transform.matrix.size().first},
      typename matrix_type::GlobalCols{transform.matrix.size().first});

  auto dualGuess = createTransformInput(transposed);
  cpppetsc::NonlinearSolver<policy_type>{
      cpppetsc::nestMatrices<policy_type>({
          {&fullStiffnessMatrix, &transposed},
          {&transform.matrix, &zeroes},
      }),
      cpppetsc::nestVectors<policy_type>({&fullForces, &dual})}
      .solve(std::move(function), std::move(jacobian),
             cpppetsc::nestVectors<policy_type>({&initialGuess, &dualGuess}));

  return initialGuess;
}

} // namespace solve
} // namespace ae108
