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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TAOSolverDivergedException.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector.h"
#include <functional>
#include <petscmat.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petsctao.h>
#include <petscvec.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class LeastSquaresSolver {
public:
  using vector_type = Vector<Policy>;
  using matrix_type = Matrix<Policy>;

  using size_type = typename Vector<Policy>::size_type;
  using value_type = typename Vector<Policy>::value_type;
  using real_type = typename Vector<Policy>::real_type;

  /**
   * @brief Initializes the nonlinear solver with the problem dimension.
   *
   * @param bufferMatrix The matrix to use as an internal buffer.
   * @param bufferVector The matrix to use as an internal buffer.
   */
  explicit LeastSquaresSolver(matrix_type bufferMatrix,
                              distributed<vector_type> bufferVector);

  /**
   * @remark The first parameter is the input, the second parameter is the
   * output.
   */
  using ObjectiveFunctor = std::function<void(const distributed<vector_type> &,
                                              distributed<vector_type> *)>;

  /**
   * @remark The first parameter is the input, the second parameter is the
   * output.
   */
  using GradientFunctor =
      std::function<void(const distributed<vector_type> &, matrix_type *)>;

  /**
   * @brief Minimizes the sum of the squared elements of the objective function.
   *
   * @param objective The function to minimize.
   * @param gradient The gradient of the objective.
   * @param guess The starting point for the optimization.
   *
   * @throw TAOSolverDivergedException if solve did not converge.
   *
   * @return The minimizer.
   */
  distributed<vector_type> solve(ObjectiveFunctor objective,
                                 GradientFunctor gradient,
                                 distributed<vector_type> guess) const;

private:
  static PetscErrorCode objectiveAdapter(Tao, Vec, Vec, void *);
  static PetscErrorCode gradientAdapter(Tao, Vec, Mat, Mat, void *);

  static Tao createTao();

  /**
   * @throw TAOSolverDivergedException if solve did not converge.
   */
  void checkConvergence() const;

  UniqueEntity<Tao> _tao;

  matrix_type _buffer_matrix;
  distributed<vector_type> _buffer_vector;
};

extern template class LeastSquaresSolver<SequentialComputePolicy>;
extern template class LeastSquaresSolver<ParallelComputePolicy>;
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cpppetsc {
template <class Policy> Tao LeastSquaresSolver<Policy>::createTao() {
  Tao tao;
  Policy::handleError(TaoCreate(Policy::communicator(), &tao));
  return tao;
}

template <class Policy>
LeastSquaresSolver<Policy>::LeastSquaresSolver(
    matrix_type bufferMatrix, distributed<vector_type> bufferVector)
    : _tao(makeUniqueEntity<Policy>(createTao())),
      _buffer_matrix(std::move(bufferMatrix)),
      _buffer_vector(std::move(bufferVector)) {
  Policy::handleError(TaoSetFromOptions(_tao.get()));
  Policy::handleError(TaoSetType(_tao.get(), TAOBRGN));
  Policy::handleError(TaoBRGNSetRegularizerWeight(_tao.get(), 0.));
}

template <class Policy>
void LeastSquaresSolver<Policy>::checkConvergence() const {
  const auto errorCode = [this]() {
    auto reason = TaoConvergedReason{};
    Policy::handleError(TaoGetConvergedReason(_tao.get(), &reason));
    return reason;
  }();

  if (errorCode < 0) {
    throw TAOSolverDivergedException{};
  }
}

template <class Policy>
distributed<typename LeastSquaresSolver<Policy>::vector_type>
LeastSquaresSolver<Policy>::solve(ObjectiveFunctor objective,
                                  GradientFunctor gradient,
                                  distributed<vector_type> guess) const {
  Policy::handleError(
      TaoSetResidualRoutine(_tao.get(), _buffer_vector.unwrap().data(),
                            &LeastSquaresSolver::objectiveAdapter, &objective));
  Policy::handleError(TaoSetJacobianResidualRoutine(
      _tao.get(), _buffer_matrix.data(), _buffer_matrix.data(),
      &LeastSquaresSolver::gradientAdapter, &gradient));
  Policy::handleError(TaoSetInitialVector(_tao.get(), guess.unwrap().data()));
  Policy::handleError(TaoSolve(_tao.get()));

  checkConvergence();

  return guess;
}

template <class Policy>
PetscErrorCode LeastSquaresSolver<Policy>::objectiveAdapter(Tao, Vec input,
                                                            Vec output,
                                                            void *context) {
  const auto functionPtr = static_cast<ObjectiveFunctor *>(context);
  distributed<vector_type> wrappedOutput(UniqueEntity<Vec>(output, [](Vec) {}));
  (*functionPtr)(distributed<vector_type>(UniqueEntity<Vec>(input, [](Vec) {})),
                 &wrappedOutput);
  return 0;
}

template <class Policy>
PetscErrorCode LeastSquaresSolver<Policy>::gradientAdapter(Tao, Vec input,
                                                           Mat output, Mat,
                                                           void *context) {
  const auto functionPtr = static_cast<GradientFunctor *>(context);
  matrix_type wrappedOutput(UniqueEntity<Mat>(output, [](Mat) {}));
  (*functionPtr)(distributed<vector_type>(UniqueEntity<Vec>(input, [](Vec) {})),
                 &wrappedOutput);
  return 0;
}

} // namespace cpppetsc
} // namespace ae108

#endif