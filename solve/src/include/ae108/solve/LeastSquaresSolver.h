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
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/clone.h"
#include <cassert>
#include <functional>

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

  /**
   * @brief Calls computeSolution with the default assembler calls (ie. calls
   * assembleForceVector or assembleStiffnessMatrix). Passes on the rest of the
   * arguments.
   *
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @para assembler Valid nonzero pointer.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(cpppetsc::distributed<vector_type> initialGuess,
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
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the local force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type>
  computeSolution(cpppetsc::distributed<vector_type> initialGuess,
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
   * E(x) -> min. More precisely, this function minimizes |F(x)|^2.
   *
   * @param initialGuess The global vector to start iterating at.
   * @param time This will be used to configure the assembler.
   * @param assembleForceVector A valid callable. It will be called to assemble
   * the distributed force vector.
   * @param assembleStiffnessMatrix A valid callable. It will be called to
   * assemble the stiffness matrix.
   */
  cpppetsc::distributed<vector_type> computeSolution(
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

template <class Assembler>
LeastSquaresSolver<Assembler>::LeastSquaresSolver(const mesh_type *mesh)
    : _mesh(mesh) {}

template <class Assembler>
cpppetsc::distributed<typename LeastSquaresSolver<Assembler>::vector_type>
LeastSquaresSolver<Assembler>::computeSolution(
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    const Assembler *const assembler) const {
  assert(assembler);
  return computeSolution(
      std::move(initialGuess), time,
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

  return computeSolution(std::move(initialGuess), time,
                         distributedForceVectorAssembler,
                         distributedStiffnessMatrixAssembler);
}

template <class Assembler>
cpppetsc::distributed<typename LeastSquaresSolver<Assembler>::vector_type>
LeastSquaresSolver<Assembler>::computeSolution(
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  auto function = [&](const cpppetsc::distributed<vector_type> &input,
                      cpppetsc::distributed<vector_type> *const out) {
    assert(out);
    out->unwrap().setZero();
    assembleForceVector(input, time, out);
  };

  auto jacobian = [&](const cpppetsc::distributed<vector_type> &input,
                      matrix_type *const out) {
    assert(out);
    out->setZero();
    assembleStiffnessMatrix(input, time, out);
  };

  const cpppetsc::LeastSquaresSolver<policy_type> solver{
      matrix_type::fromMesh(*_mesh), vector_type::fromGlobalMesh(*_mesh)};

  return solver.solve(std::move(function), std::move(jacobian),
                      std::move(initialGuess));
}

} // namespace solve
} // namespace ae108

#endif