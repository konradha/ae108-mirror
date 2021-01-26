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
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/clone.h"
#include "ae108/solve/NonlinearSolver.h"
#include "ae108/solve/dynamics/DynamicState.h"
#include "ae108/solve/dynamics/NewmarkParameters.h"
#include <cassert>
#include <functional>

namespace ae108 {
namespace solve {

template <class Assembler, class NonlinearSolver = NonlinearSolver<Assembler>>
class DynamicSolver {
public:
  using mesh_type = typename assembly::MeshTypeTrait<Assembler>::type;
  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;

  using nonlinearsolver_type = NonlinearSolver;
  using state_type = dynamics::DynamicState<cpppetsc::distributed<vector_type>>;

  /**
   * @param mesh A valid pointer to a mesh_type instance.
   * @param solver A valid pointer to a nonlinearsolver_type instance.
   */
  explicit DynamicSolver(const mesh_type *mesh,
                         const nonlinearsolver_type *solver,
                         const dynamics::NewmarkParameters &newmark);

  using BoundaryConditionContainer =
      typename nonlinearsolver_type::BoundaryConditionContainer;

  /**
   * @brief Calls another overload of computeSolution passing on all arguments
   * except assembler. The local force vector and the local stiffness matrix are
   * computed using assembler's assembleForceVector and assembleStiffnessMatrix
   * methods.
   *
   * @param boundaryConditions The local essential boundary conditions to apply.
   * @param state The state to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param timestep The current timestep.
   * @param mass The mass matrix.
   * @param damping The damping matrix.
   * @param assembler Valid nonzero pointer.
   */
  state_type
  computeSolution(const BoundaryConditionContainer &boundaryConditions,
                  state_type state, const double time, const double timestep,
                  const matrix_type &mass, const matrix_type &damping,
                  const Assembler *const assembler) const;

  using LocalForceVectorAssembler =
      std::function<void(const cpppetsc::local<vector_type> &, double,
                         cpppetsc::local<vector_type> *)>;

  using LocalStiffnessMatrixAssembler = std::function<void(
      const cpppetsc::local<vector_type> &, double, matrix_type *)>;

  /**
   * @brief Calls another overload of computeSolution passing on all arguments
   * except the assembling functions. The global force vector and the global
   * stiffness matrix are computed using the local contributions.
   *
   * @param boundaryConditions The local essential boundary conditions to apply.
   * @param state The state to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param timestep The current timestep.
   * @param mass The mass matrix.
   * @param damping The damping matrix.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the local force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the local stiffness matrix.
   */
  state_type
  computeSolution(const BoundaryConditionContainer &boundaryConditions,
                  state_type state, const double time, const double timestep,
                  const matrix_type &mass, const matrix_type &damping,
                  LocalForceVectorAssembler assembleForceVector,
                  LocalStiffnessMatrixAssembler assembleStiffnessMatrix) const;

  using DistributedForceVectorAssembler =
      std::function<void(const cpppetsc::distributed<vector_type> &, double,
                         cpppetsc::distributed<vector_type> *)>;

  using DistributedStiffnessMatrixAssembler = std::function<void(
      const cpppetsc::distributed<vector_type> &, double, matrix_type *)>;

  /**
   * @brief Computes the solution using the provided nonlinear solver type.
   *
   * @param boundaryConditions The local essential boundary conditions to apply.
   * @param state The state to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param timestep The current timestep.
   * @param mass The mass matrix.
   * @param damping The damping matrix.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the global force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the global stiffness matrix.
   */
  state_type computeSolution(
      const BoundaryConditionContainer &boundaryConditions, state_type state,
      const double time, const double timestep, const matrix_type &mass,
      const matrix_type &damping,
      DistributedForceVectorAssembler assembleForceVector,
      DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const;

private:
  const mesh_type *_mesh;
  const nonlinearsolver_type *_solver;
  dynamics::NewmarkParameters _newmark;

  /*
   * @brief Computes the vector that only depends on the current state:
   * M * (
   *         1 / (beta * dt^2) * U +
   *         1 / (beta * dt) * U' +
   *         1 / (2 * beta -1) * U''
   *     ) +
   * C * (
   *         gamma / (beta * dt) * U +
   *         (gamma/beta - 1) * U' +
   *         dt * (gamma / (2 * beta) - 1) * U''
   *     )
   */
  cpppetsc::distributed<vector_type>
  computeRhsVector(const state_type &state, const double timestep,
                   const matrix_type &mass, const matrix_type &damping) const;

  /*
   * @brief Computes the matrix that only depends on the current state:
   * 1 / (beta * dt^2) * M + gamma / (beta * dt) * C
   */
  matrix_type computeLhsMatrix(const double timestep, const matrix_type &mass,
                               const matrix_type &damping) const;
};
} // namespace solve
} // namespace ae108

/********************************************************************
 *  implementations of methods
 *******************************************************************/

namespace ae108 {
namespace solve {

template <class Assembler, class NonlinearSolver>
DynamicSolver<Assembler, NonlinearSolver>::DynamicSolver(
    const mesh_type *mesh, const nonlinearsolver_type *const solver,
    const dynamics::NewmarkParameters &newmark)
    : _mesh(mesh), _solver(solver), _newmark(newmark) {
  assert(_mesh);
  assert(_solver);
}

template <class Assembler, class NonlinearSolver>
typename DynamicSolver<Assembler, NonlinearSolver>::state_type
DynamicSolver<Assembler, NonlinearSolver>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions, state_type state,
    const double time, const double timestep, const matrix_type &mass,
    const matrix_type &damping, const Assembler *const assembler) const {
  assert(assembler);
  return computeSolution(
      boundaryConditions, std::move(state), time, timestep, mass, damping,
      [assembler](const cpppetsc::local<vector_type> &displacements,
                  const double time,
                  cpppetsc::local<vector_type> *const output) {
        assembler->assembleForceVector(displacements, time, output);
      },
      [assembler](const cpppetsc::local<vector_type> &displacements,
                  const double time, matrix_type *const output) {
        assembler->assembleStiffnessMatrix(displacements, time, output);
      });
}

template <class Assembler, class NonlinearSolver>
typename DynamicSolver<Assembler, NonlinearSolver>::state_type
DynamicSolver<Assembler, NonlinearSolver>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions, state_type state,
    const double time, const double timestep, const matrix_type &mass,
    const matrix_type &damping, LocalForceVectorAssembler assembleForceVector,
    LocalStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  const auto &mesh = *_mesh;
  auto localDisplacements = vector_type::fromLocalMesh(mesh);
  auto localForces = vector_type::fromLocalMesh(mesh);

  return computeSolution(
      boundaryConditions, std::move(state), time, timestep, mass, damping,
      [&localDisplacements, &localForces, &mesh, &assembleForceVector,
       &assembleStiffnessMatrix](
          const cpppetsc::distributed<vector_type> &displacements,
          const double time, cpppetsc::distributed<vector_type> *const output) {
        mesh.copyToLocalVector(displacements, &localDisplacements);
        localForces.unwrap().setZero();
        assembleForceVector(localDisplacements, time, &localForces);
        mesh.addToGlobalVector(localForces, output);
      },
      [&localDisplacements, &mesh, &assembleForceVector,
       assembleStiffnessMatrix](
          const cpppetsc::distributed<vector_type> &displacements,
          const double time, matrix_type *const output) {
        mesh.copyToLocalVector(displacements, &localDisplacements);
        assembleStiffnessMatrix(localDisplacements, time, output);
        output->finalize();
      });
}

template <class Assembler, class NonlinearSolver>
typename DynamicSolver<Assembler, NonlinearSolver>::state_type
DynamicSolver<Assembler, NonlinearSolver>::computeSolution(
    const BoundaryConditionContainer &boundaryConditions, state_type state,
    const double time, const double timestep, const matrix_type &mass,
    const matrix_type &damping,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  const auto lhsMatrix = computeLhsMatrix(timestep, mass, damping);
  const auto rhsVector = computeRhsVector(state, timestep, mass, damping);

  state_type newState = {clone(state.displacements), clone(state.velocities),
                         clone(state.accelerations)};

  const auto &mesh = *_mesh;
  auto localDisplacements = vector_type::fromLocalMesh(mesh);
  auto localForces = vector_type::fromLocalMesh(mesh);

  using ForceVectorAssembler =
      typename nonlinearsolver_type::DistributedForceVectorAssembler;
  using StiffnessMatrixAssembler =
      typename nonlinearsolver_type::DistributedStiffnessMatrixAssembler;

  newState.displacements = _solver->computeSolution(
      boundaryConditions, std::move(newState.displacements), time,
      ForceVectorAssembler(
          [&](const cpppetsc::distributed<vector_type> &input,
              const double time,
              cpppetsc::distributed<vector_type> *const output) {
            assert(output);

            assembleForceVector(input, time, output);
            output->unwrap().addAx(lhsMatrix, input);
            output->unwrap().timesAlphaPlusBetaX(1., -1., rhsVector);
          }),
      StiffnessMatrixAssembler(
          [&](const cpppetsc::distributed<vector_type> &input,
              const double time, matrix_type *const output) {
            assert(output);

            assembleStiffnessMatrix(input, time, output);
            output->addAlphaX(1., lhsMatrix);
          }));

  // newState.displacements now contains final result

  state.displacements.unwrap().timesAlphaPlusBetaXPlusGammaY(
      -1., 1., newState.displacements, -1. * timestep, state.velocities);
  const auto factor_1 = 1. / (_newmark.beta * timestep * timestep);
  const auto factor_2 = 1. - 1. / (2. * _newmark.beta);
  newState.accelerations.unwrap().timesAlphaPlusBetaX(factor_2, factor_1,
                                                      state.displacements);

  // newState.accelerations now contains final result
  // state.displacements now contains invalid value

  newState.velocities.unwrap().timesAlphaPlusBetaXPlusGammaY(
      1., timestep * _newmark.gamma, newState.accelerations,
      timestep * (1. - _newmark.gamma), state.accelerations);

  // newState.velocities now contains final result

  return newState;
}

template <class Assembler, class NonlinearSolver>
cpppetsc::distributed<
    typename DynamicSolver<Assembler, NonlinearSolver>::vector_type>
DynamicSolver<Assembler, NonlinearSolver>::computeRhsVector(
    const state_type &state, const double timestep, const matrix_type &mass,
    const matrix_type &damping) const {
  auto result = vector_type::fromLayoutOf(state.displacements);

  // add M * (
  //         1 / (beta * dt^2) * U +
  //         1 / (beta * dt) * U' +
  //         1 / (2 * beta -1) * U''
  //         )
  result.addAx(mass, [&]() {
    auto vec = clone(state.displacements);
    vec.unwrap().timesAlphaPlusBetaXPlusGammaY(
        1. / _newmark.beta / timestep / timestep, 1. / _newmark.beta / timestep,
        state.velocities, (1. / 2. / _newmark.beta - 1.), state.accelerations);
    return vec;
  }());

  // add C * (
  //         gamma / (beta * dt) * U +
  //         (gamma/beta - 1) * U' +
  //         dt * (gamma / (2 * beta) - 1) * U''
  //         )
  result.addAx(damping, [&]() {
    auto vec = clone(state.displacements);
    vec.unwrap().timesAlphaPlusBetaXPlusGammaY(
        _newmark.gamma / _newmark.beta / timestep,
        (_newmark.gamma / _newmark.beta - 1.), state.velocities,
        timestep * (_newmark.gamma / 2. / _newmark.beta - 1.),
        state.accelerations);
    return vec;
  }());

  return cpppetsc::tag<cpppetsc::DistributedTag>(std::move(result));
}

template <class Assembler, class NonlinearSolver>
typename DynamicSolver<Assembler, NonlinearSolver>::matrix_type
DynamicSolver<Assembler, NonlinearSolver>::computeLhsMatrix(
    const double timestep, const matrix_type &mass,
    const matrix_type &damping) const {
  auto matrix = matrix_type::clone(mass);

  // compute 1 / (beta * dt^2) * M
  matrix.scale(1. / _newmark.beta / timestep / timestep);

  // add gamma / (beta * dt) * C
  matrix.addAlphaX(_newmark.gamma / _newmark.beta / timestep, damping);

  return matrix;
}

} // namespace solve
} // namespace ae108
