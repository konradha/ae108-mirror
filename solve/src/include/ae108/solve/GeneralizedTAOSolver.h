// © 2020, 2022 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/cpppetsc/MeshBoundaryCondition.h"
#include "ae108/cpppetsc/TAOSolver.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include <cassert>
#include <functional>
#include <mpi.h>
#include <petscsys.h>
#include <vector>

namespace ae108 {
namespace solve {

template <class Assembler> class GeneralizedTAOSolver {
public:
  using policy_type = typename assembly::PolicyTypeTrait<Assembler>::type;
  using mesh_type = typename assembly::MeshTypeTrait<Assembler>::type;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;
  using real_type = typename mesh_type::real_type;
  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;

  /**
   * @param mesh A valid nonzero pointer to a mesh_type instance.
   */
  explicit GeneralizedTAOSolver(const mesh_type *mesh);

  struct BoundaryConditionContainer {
    /**
     * @brief The dimension of the output of `residual`.
     */
    size_type size;

    /**
     * @brief The solver attempts to reach a residual of zero.
     */
    typename cpppetsc::TAOSolver<policy_type>::GradientFunctor residual;

    /**
     * @brief The jacobian of the residual.
     */
    typename cpppetsc::TAOSolver<policy_type>::HessianFunctor jacobian;
  };

  /**
   * @brief Calls computeSolution with the default assembler calls (ie. calls
   * assembleEnergy, assembleForceVector or assembleStiffnessMatrix). Passes on
   * the rest of the arguments.
   *
   * @param boundaryConditions The nonlinear boundary conditions to apply.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembler Valid nonzero pointer.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(const BoundaryConditionContainer &boundaryConditions,
                  cpppetsc::distributed<vector_type> initialGuess,
                  const double time, const Assembler *const assembler) const;

  using LocalEnergyAssembler = std::function<void(
      const cpppetsc::local<vector_type> &, double, real_type *)>;

  using LocalForceVectorAssembler =
      std::function<void(const cpppetsc::local<vector_type> &, double,
                         cpppetsc::local<vector_type> *)>;

  using LocalStiffnessMatrixAssembler = std::function<void(
      const cpppetsc::local<vector_type> &, double, matrix_type *)>;

  /**
   * @brief Calls computeSolution wrapping up the local assembler calls into
   * distributed assembler calls. Passes on the rest of the arguments.
   *
   * @param boundaryConditions The nonlinear boundary conditions to apply.
   * Note that empty functions (residual, jacobian) are interpreted as functions
   * that do not change the output.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleEnergy A valid callable. It will be called to assemble
   * the local energy.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the local force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(const BoundaryConditionContainer &boundaryConditions,
                  cpppetsc::distributed<vector_type> initialGuess,
                  const double time, LocalEnergyAssembler assembleEnergy,
                  LocalForceVectorAssembler assembleForceVector,
                  LocalStiffnessMatrixAssembler assemblerStiffnessMatrix) const;

  using DistributedEnergyAssembler = std::function<void(
      const cpppetsc::distributed<vector_type> &, double, real_type *)>;

  using DistributedForceVectorAssembler =
      std::function<void(const cpppetsc::distributed<vector_type> &, double,
                         cpppetsc::distributed<vector_type> *)>;

  using DistributedStiffnessMatrixAssembler = std::function<void(
      const cpppetsc::distributed<vector_type> &, double, matrix_type *)>;

  /**
   * @brief The goal of this function is to solve E(x) = min under the condition
   * that the residual is equal to zero.
   *
   * @remark Does not update internal variables.
   *
   * @param boundaryConditions The nonlinear boundary conditions to apply.
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleEnergy A valid callable. It will be called to assemble
   * the distributed energy.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the distributed force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type> computeSolution(
      const BoundaryConditionContainer &boundaryConditions,
      cpppetsc::distributed<vector_type> initialGuess, const double time,
      DistributedEnergyAssembler assembleEnergy,
      DistributedForceVectorAssembler assembleForceVector,
      DistributedStiffnessMatrixAssembler assemblerStiffnessMatrix) const;

private:
  const mesh_type *_mesh;
}; // namespace solve
} // namespace solve
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace solve {
template <class Assembler>
GeneralizedTAOSolver<Assembler>::GeneralizedTAOSolver(const mesh_type *mesh)
    : _mesh(mesh) {}

template <class Assembler>
cpppetsc::distributed<typename GeneralizedTAOSolver<Assembler>::vector_type>
GeneralizedTAOSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    const Assembler *const assembler) const {
  assert(assembler);
  return computeSolution(
      boundaryConditions, std::move(initialGuess), time,
      LocalEnergyAssembler(
          [assembler](const cpppetsc::local<vector_type> &localDisplacements,
                      const double time, real_type *const localEnergy) {
            assembler->assembleEnergy(localDisplacements, time, localEnergy);
          }),
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
cpppetsc::distributed<typename GeneralizedTAOSolver<Assembler>::vector_type>
GeneralizedTAOSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    LocalEnergyAssembler assembleEnergy,
    LocalForceVectorAssembler assembleForceVector,
    LocalStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleEnergy);
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  const auto &mesh = *_mesh;

  auto localDisplacements = vector_type::fromLocalMesh(mesh);
  auto localForces = vector_type::fromLocalMesh(mesh);

  const auto distributedEnergyAssembler =
      [&assembleEnergy, &mesh,
       &localDisplacements](const cpppetsc::distributed<vector_type> &input,
                            const double time, real_type *const output) {
        mesh.copyToLocalVector(input, &localDisplacements);

        *output = real_type{};
        assembleEnergy(localDisplacements, time, output);

        policy_type::handleError(MPI_Allreduce(MPI_IN_PLACE, output, 1,
                                               MPIU_REAL, MPIU_SUM,
                                               policy_type::communicator()));
      };

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
                         distributedEnergyAssembler,
                         distributedForceVectorAssembler,
                         distributedStiffnessMatrixAssembler);
}

template <class Assembler>
cpppetsc::distributed<typename GeneralizedTAOSolver<Assembler>::vector_type>
GeneralizedTAOSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    DistributedEnergyAssembler assembleEnergy,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleEnergy);
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  using Solver = cpppetsc::TAOSolver<policy_type>;
  auto matrix = matrix_type::fromMesh(*_mesh);

  const auto global = matrix.size();
  const auto local = matrix.localSize();

  Solver solver{std::move(matrix), Solver::Type::pdipm};

  solver.setConstraints(typename Solver::EqualityConstraints{
      boundaryConditions.residual
          ? boundaryConditions.residual
          : [](const cpppetsc::distributed<vector_type> &,
               cpppetsc::distributed<vector_type> *) {},
      boundaryConditions.jacobian
          ? boundaryConditions.jacobian
          : [](const cpppetsc::distributed<vector_type> &, matrix_type *) {},
      {
          cpppetsc::distributed<vector_type>(boundaryConditions.size),
          matrix_type(typename matrix_type::LocalRows(PETSC_DECIDE),
                      typename matrix_type::LocalCols(local.second),
                      typename matrix_type::GlobalRows(boundaryConditions.size),
                      typename matrix_type::GlobalCols(global.second)),
      },
  });

  return solver.solve(
      [time, &assembleEnergy](const cpppetsc::distributed<vector_type> &input,
                              double *const output) {
        *output = real_type{};
        assembleEnergy(input, time, output);
      },
      [time,
       &assembleForceVector](const cpppetsc::distributed<vector_type> &input,
                             cpppetsc::distributed<vector_type> *const output) {
        output->unwrap().setZero();
        assembleForceVector(input, time, output);
      },
      [time, &assembleStiffnessMatrix](
          const cpppetsc::distributed<vector_type> &input,
          matrix_type *const output) {
        output->setZero();
        assembleStiffnessMatrix(input, time, output);
      },
      std::move(initialGuess));
}

} // namespace solve
} // namespace ae108

#endif