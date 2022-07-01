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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/NonlinearSolverDivergedException.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector.h"
#include <functional>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscvec.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class NonlinearSolver {
public:
  using vector_type = Vector<Policy>;
  using matrix_type = Matrix<Policy>;

  using size_type = typename Matrix<Policy>::size_type;

  /**
   * @brief Initializes the nonlinear solver with the problem dimension.
   *
   * @param dimension The ambient dimension, e.g. the size of the vector that
   * the function returns.
   */
  explicit NonlinearSolver(const size_type dimension);

  /**
   * @brief Initializes the nonlinear solver with the given buffer entities of
   * correct dimension.
   *
   * @remark If `buffer_matrix` is of type MATNEST or MATSHELL, then no
   * preconditioner is applied by default.
   */
  explicit NonlinearSolver(matrix_type buffer_matrix,
                           distributed<vector_type> buffer_vector);

  /**
   * @remark The first parameter is the input, the second parameter is the
   * output.
   */
  using FunctionFunctor = std::function<void(const distributed<vector_type> &,
                                             distributed<vector_type> *)>;

  /**
   * @brief The first parameter is the input, the second parameter is the
   * output.
   */
  using JacobianFunctor =
      std::function<void(const distributed<vector_type> &, matrix_type *)>;

  /**
   * @brief Solves the nonlinear equation f(x) == 0 for x.
   *
   * @param function A function object that computes f at a given point.
   * @param jacobian A function object that computes the Jacobian of f at
   * given point.
   *
   * @throw NonlinearSolverDivergedException if the solve diverged.
   *
   * @return The solution x of the nonlinear equation.
   */
  distributed<vector_type> solve(FunctionFunctor function,
                                 JacobianFunctor jacobian,
                                 distributed<vector_type> guess) const;

private:
  static PetscErrorCode functionAdapter(SNES, Vec, Vec, void *);
  static PetscErrorCode jacobianAdapter(SNES, Vec, Mat, Mat, void *);

  static SNES createSNES();

  UniqueEntity<SNES> _snes;

  distributed<vector_type> _buffer_vector;
  matrix_type _buffer_matrix;
};

extern template class NonlinearSolver<SequentialComputePolicy>;
extern template class NonlinearSolver<ParallelComputePolicy>;
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cpppetsc {

template <class Policy> SNES NonlinearSolver<Policy>::createSNES() {
  SNES snes;
  Policy::handleError(SNESCreate(Policy::communicator(), &snes));
  return snes;
}

template <class Policy>
NonlinearSolver<Policy>::NonlinearSolver(const size_type dimension)
    : _snes(makeUniqueEntity<Policy>(createSNES())), _buffer_vector(dimension),
      _buffer_matrix(dimension, dimension) {
  Policy::handleError(SNESSetFromOptions(_snes.get()));
}

template <class Policy>
NonlinearSolver<Policy>::NonlinearSolver(matrix_type buffer_matrix,
                                         distributed<vector_type> buffer_vector)
    : _snes(makeUniqueEntity<Policy>(createSNES())),
      _buffer_vector(std::move(buffer_vector)),
      _buffer_matrix(std::move(buffer_matrix)) {
  const auto deactivatePreconditioning = [&]() {
    auto ksp = KSP{};
    Policy::handleError(SNESGetKSP(_snes.get(), &ksp));
    auto pc = PC{};
    Policy::handleError(KSPGetPC(ksp, &pc));
    Policy::handleError(PCSetType(pc, PCNONE));
  };

  auto type = MatType{};
  Policy::handleError(MatGetType(_buffer_matrix.data(), &type));

  if (std::string(type) == std::string(MATNEST) ||
      std::string(type) == std::string(MATSHELL)) {
    deactivatePreconditioning();
  }
  Policy::handleError(SNESSetFromOptions(_snes.get()));
}

template <class Policy>
distributed<typename NonlinearSolver<Policy>::vector_type>
NonlinearSolver<Policy>::solve(FunctionFunctor function,
                               JacobianFunctor jacobian,
                               distributed<vector_type> guess) const {
  Policy::handleError(
      SNESSetFunction(_snes.get(), _buffer_vector.unwrap().data(),
                      &NonlinearSolver::functionAdapter, &function));
  Policy::handleError(
      SNESSetJacobian(_snes.get(), _buffer_matrix.data(), _buffer_matrix.data(),
                      &NonlinearSolver::jacobianAdapter, &jacobian));
  Policy::handleError(SNESSolve(_snes.get(), nullptr, guess.unwrap().data()));

  const auto errorCode = [this]() {
    auto code = SNESConvergedReason{};
    Policy::handleError(SNESGetConvergedReason(_snes.get(), &code));
    return code;
  }();

  if (errorCode < 0) {
    throw NonlinearSolverDivergedException{};
  }

  return guess;
}

template <class Policy>
PetscErrorCode NonlinearSolver<Policy>::functionAdapter(SNES, Vec input,
                                                        Vec output,
                                                        void *context) {
  const auto functionPtr = static_cast<FunctionFunctor *>(context);
  distributed<vector_type> wrappedOutput(UniqueEntity<Vec>(output, [](Vec) {}));
  (*functionPtr)(distributed<vector_type>(UniqueEntity<Vec>(input, [](Vec) {})),
                 &wrappedOutput);
  return 0;
}

template <class Policy>
PetscErrorCode NonlinearSolver<Policy>::jacobianAdapter(SNES, Vec input, Mat,
                                                        Mat output,
                                                        void *context) {
  const auto jacobianPtr = static_cast<JacobianFunctor *>(context);
  matrix_type wrappedOutput(UniqueEntity<Mat>(output, [](Mat) {}));
  (*jacobianPtr)(distributed<vector_type>(UniqueEntity<Vec>(input, [](Vec) {})),
                 &wrappedOutput);
  return 0;
}
} // namespace cpppetsc
} // namespace ae108
