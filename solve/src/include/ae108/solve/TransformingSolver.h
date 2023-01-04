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
#include "ae108/cpppetsc/clone.h"
#include "ae108/solve/AffineTransform.h"
#include <cassert>
#include <functional>
#include <vector>

namespace ae108 {
namespace solve {

template <class Assembler> class TransformingSolver {
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
  explicit TransformingSolver(const mesh_type *mesh);

  using BoundaryConditionContainer = AffineTransform<policy_type>;

  /**
   * @brief Calls computeSolution with the default assembler calls (ie. calls
   * assembleForceVector or assembleStiffnessMatrix). Passes on the rest of the
   * arguments.
   *
   * @param transform The affine transform to apply.
   * @param initialGuess The global (untransformed) vector to start iterating
   * at.
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
   * @param transform The affine transform to apply.
   * @param initialGuess The global (untransformed) vector to start iterating
   * at.
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
   * @brief The goal of this function is to find the minimizer y of
   * E(P * y + c) -> min, where (P, c) is the provided affine transformation
   * `transform`. The transform can be used to specify boundary conditions by
   * fixing degrees of freedom in the output of the transform.
   *
   * More precisely, this function solves the equation E' = 0, which is
   * equivalent to P^t * F(P * y + c) = 0, where F(P * y + c) are the forces at
   * P * y + c.
   *
   * @param transform The affine transform to apply.
   * @param initialGuess The global (untransformed) vector to start iterating
   * at.
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
TransformingSolver<Assembler>::TransformingSolver(const mesh_type *mesh)
    : _mesh(mesh) {}

template <class Assembler>
cpppetsc::distributed<typename TransformingSolver<Assembler>::vector_type>
TransformingSolver<Assembler>::computeSolution(
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
cpppetsc::distributed<typename TransformingSolver<Assembler>::vector_type>
TransformingSolver<Assembler>::computeSolution(
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
cpppetsc::distributed<typename TransformingSolver<Assembler>::vector_type>
TransformingSolver<Assembler>::computeSolution(
    const BoundaryConditionContainer &transform,
    cpppetsc::distributed<vector_type> initialGuess, const double time,
    DistributedForceVectorAssembler assembleForceVector,
    DistributedStiffnessMatrixAssembler assembleStiffnessMatrix) const {
  assert(assembleForceVector);
  assert(assembleStiffnessMatrix);

  auto fullInput = clone(transform.vector);

  auto fullForces = vector_type::fromGlobalMesh(*_mesh);
  auto function = [&](const cpppetsc::distributed<vector_type> &input,
                      cpppetsc::distributed<vector_type> *const output) {
    apply(transform, input, &fullInput);

    fullForces.unwrap().setZero();
    assembleForceVector(fullInput, time, &fullForces);
    policy_type::handleError(MatMultTranspose(transform.matrix.data(),
                                              fullForces.unwrap().data(),
                                              output->unwrap().data()));
  };

  auto fullStiffnessMatrix = matrix_type::fromMesh(*_mesh);
  auto jacobian = [&](const cpppetsc::distributed<vector_type> &input,
                      matrix_type *const output) {
    apply(transform, input, &fullInput);

    fullStiffnessMatrix.setZero();
    assembleStiffnessMatrix(fullInput, time, &fullStiffnessMatrix);
    auto outputData = output->data();
    policy_type::handleError(MatPtAP(fullStiffnessMatrix.data(),
                                     transform.matrix.data(), MAT_REUSE_MATRIX,
                                     PETSC_DEFAULT, &outputData));
    assert(outputData == output->data() && "PETSc does not change the pointer");
  };

  const cpppetsc::NonlinearSolver<policy_type> solver{
      matrix_type::fromPtAP(transform.matrix, fullStiffnessMatrix),
      createTransformInput(transform.matrix)};
  return solver.solve(std::move(function), std::move(jacobian),
                      std::move(initialGuess));
}

} // namespace solve
} // namespace ae108
