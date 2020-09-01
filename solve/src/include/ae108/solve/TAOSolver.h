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

template <class Assembler> class TAOSolver {
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
  explicit TAOSolver(const mesh_type *mesh);

  using BoundaryConditionContainer =
      std::vector<cpppetsc::MeshBoundaryCondition<mesh_type>>;

  /**
   * @brief Calls computeSolution with the default assembler calls (ie. calls
   * assembleEnergy, assembleForceVector or assembleStiffnessMatrix). Passes on
   * the rest of the arguments.
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

  using LocalEnergyAssembler = std::function<void(
      const cpppetsc::local<vector_type> &, double, value_type *)>;

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
      const cpppetsc::distributed<vector_type> &, double, value_type *)>;

  using DistributedForceVectorAssembler =
      std::function<void(const cpppetsc::distributed<vector_type> &, double,
                         cpppetsc::distributed<vector_type> *)>;

  using DistributedStiffnessMatrixAssembler = std::function<void(
      const cpppetsc::distributed<vector_type> &, double, matrix_type *)>;

  /**
   * @brief The goal of this function is to solve E(x) = min under the condition
   * that all the essential boundary conditions x_i = a_i are met.
   *
   * @remark Does not update internal variables.
   *
   * @param boundaryConditions The local essential boundary conditions to apply.
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
};
} // namespace solve
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace solve {
template <class Assembler>
TAOSolver<Assembler>::TAOSolver(const mesh_type *mesh) : _mesh(mesh) {}

template <class Assembler>
cpppetsc::distributed<typename TAOSolver<Assembler>::vector_type>
TAOSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    const Assembler *const assembler) const {
  assert(assembler);
  return computeSolution(
      boundaryConditions, std::move(initialGuess), time,
      LocalEnergyAssembler(
          [assembler](const cpppetsc::local<vector_type> &localDisplacements,
                      const double time, value_type *const localEnergy) {
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
cpppetsc::distributed<typename TAOSolver<Assembler>::vector_type>
TAOSolver<Assembler>::computeSolution(
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
                            const double time, value_type *const output) {
        mesh.copyToLocalVector(input, &localDisplacements);

        *output = value_type{};
        assembleEnergy(localDisplacements, time, output);

        policy_type::handleError(MPI_Allreduce(MPI_IN_PLACE, output, 1,
                                               MPIU_SCALAR, MPIU_SUM,
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
cpppetsc::distributed<typename TAOSolver<Assembler>::vector_type>
TAOSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    DistributedEnergyAssembler assembleEnergy,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleEnergy);
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  using Solver = cpppetsc::TAOSolver<policy_type>;
  auto lowerBound = vector_type::fromLocalMesh(*_mesh);
  lowerBound.unwrap().fill(Solver::no_lower_bound);

  auto upperBound = vector_type::fromLocalMesh(*_mesh);
  upperBound.unwrap().fill(Solver::no_upper_bound);

  for (const auto &bc : boundaryConditions) {
    const auto applyBound = [&bc](cpppetsc::local<vector_type> *const bound) {
      assert(bound);

      std::vector<value_type> data;
      bc.node.copyVertexData(*bound, &data);
      data.at(bc.dof) = bc.value;
      bc.node.setVertexData(data, bound);
    };

    applyBound(&lowerBound);
    applyBound(&upperBound);
  }

  Solver solver{matrix_type::fromMesh(*_mesh), Solver::Type::bnls};

  const auto distributeBounds =
      [this](const cpppetsc::local<vector_type> &localBounds) {
        auto bounds = vector_type::fromGlobalMesh(*_mesh);
        _mesh->copyToGlobalVector(localBounds, &bounds);
        return bounds;
      };

  solver.setBounds(distributeBounds(lowerBound), distributeBounds(upperBound));

  return solver.solve(
      [time, &assembleEnergy](const cpppetsc::distributed<vector_type> &input,
                              double *const output) {
        *output = value_type{};
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
