// © 2020, 2022 ETH Zurich, Mechanics and Materials Lab
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

#ifndef AE108_PETSC_COMPLEX

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TAOSolverDivergedException.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector.h"
#include <functional>
#include <optional>
#include <petscmat.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petsctao.h>
#include <petscvec.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class TAOSolver {
public:
  using vector_type = Vector<Policy>;
  using matrix_type = Matrix<Policy>;

  using size_type = typename Vector<Policy>::size_type;
  using value_type = typename Vector<Policy>::value_type;
  using real_type = typename Vector<Policy>::real_type;

  /**
   * @brief Initializes the nonlinear solver with the problem dimension.
   *
   * @param The matrix to use as an internal buffer.
   */
  explicit TAOSolver(matrix_type bufferMatrix);

  enum class Type {
    lmvm,
    nls,
    ntr,
    ntl,
    cg,
    tron,
    owlqn,
    bmrm,
    blmvm,
    bqnls,
    bncg,
    bnls,
    bntr,
    bqnkls,
    bqnktr,
    bqnktl,
    bqpip,
    gpcg,
    nm,
    pounders,
    brgn,
    lcl,
    ssils,
    ssfls,
    asils,
    asfls,
    ipm,
    pdipm,
    shell,
    admm,
    almm,
  };

  /**
   * @brief Initializes the nonlinear solver with the problem dimension.
   *
   * @param The matrix to use as an internal buffer.
   * @param defaultType The default type of solver used. This type will not be
   * used if e.g. a command line option specifying the type is present (see
   * PETSc documentation).
   */
  explicit TAOSolver(matrix_type bufferMatrix, const Type defaultType);

  /**
   * @brief Set the type of the optimizer.
   *
   * @param type The type of the optimizer. See the definition of Type and the
   * PETSc documentation for more information about the possible choices.
   */
  void setType(const Type type);

  /**
   * @remark The first parameter is the input, the second parameter is the
   * output.
   */
  using ObjectiveFunctor =
      std::function<void(const distributed<vector_type> &, real_type *)>;

  /**
   * @remark The first parameter is the input, the second parameter is the
   * output.
   */
  using GradientFunctor = std::function<void(const distributed<vector_type> &,
                                             distributed<vector_type> *)>;

  /**
   * @remark The first parameter is the input, the second parameter is the
   * output.
   */
  using HessianFunctor =
      std::function<void(const distributed<vector_type> &, matrix_type *)>;

  static constexpr value_type no_upper_bound = PETSC_INFINITY;
  static constexpr value_type no_lower_bound = PETSC_NINFINITY;

  /**
   * @brief Define upper and lower bounds for the variables.
   *
   * @param lowerBounds A vector of constraints for every variable. Use
   * no_lower_bound to set no constraints for a variable.
   * @param upperBounds A vector of constraints for every variable. Use
   * no_upper_bound to set no constraints for a variable.
   */
  void setBounds(const distributed<vector_type> &lowerBounds,
                 const distributed<vector_type> &upperBounds);

  struct EqualityConstraints {
    /**
     * @brief The solver attempts to reach a residual of zero.
     */
    GradientFunctor residual;

    /**
     * @brief The jacobian of the residual.
     */
    HessianFunctor jacobian;

    /**
     * @brief Used to store residuals and jacobians internally.
     */
    std::pair<distributed<vector_type>, matrix_type> buffer;
  };

  /**
   * @brief Set the equality constraints of the problem.
   * @note Not all solvers consider these constraints.
   */
  void setConstraints(EqualityConstraints constraints);

  /**
   * @brief Minimizes the objective function using optimization.
   *
   * @param objective The function to minimize.
   * @param gradient The gradient of the objective.
   * @param hessian The hessian of the objective.
   * @param guess The starting point for the optimization.
   *
   * @throw TAOSolverDivergedException if solve did not converge.
   *
   * @return The minimizer.
   */
  distributed<vector_type> solve(ObjectiveFunctor objective,
                                 GradientFunctor gradient,
                                 HessianFunctor hessian,
                                 distributed<vector_type> guess) const;

private:
  static PetscErrorCode objectiveAdapter(Tao, Vec, real_type *, void *);
  static PetscErrorCode gradientAdapter(Tao, Vec, Vec, void *);
  static PetscErrorCode hessianAdapter(Tao, Vec, Mat, Mat, void *);

  static Tao createTao();

  /**
   * @throw TAOSolverDivergedException if solve did not converge.
   */
  void checkConvergence() const;

  UniqueEntity<Tao> _tao;

  matrix_type _buffer_matrix;

  std::optional<EqualityConstraints> _constraints;
};

template <class Policy>
constexpr
    typename TAOSolver<Policy>::value_type TAOSolver<Policy>::no_upper_bound;
template <class Policy>
constexpr
    typename TAOSolver<Policy>::value_type TAOSolver<Policy>::no_lower_bound;

extern template class TAOSolver<SequentialComputePolicy>;
extern template class TAOSolver<ParallelComputePolicy>;
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cpppetsc {
template <class Policy> Tao TAOSolver<Policy>::createTao() {
  Tao tao;
  Policy::handleError(TaoCreate(Policy::communicator(), &tao));
  return tao;
}

template <class Policy>
TAOSolver<Policy>::TAOSolver(matrix_type bufferMatrix)
    : _tao(makeUniqueEntity<Policy>(createTao())),
      _buffer_matrix(std::move(bufferMatrix)) {
  Policy::handleError(TaoSetFromOptions(_tao.get()));
}

template <class Policy>
TAOSolver<Policy>::TAOSolver(matrix_type bufferMatrix, const Type type)
    : _tao(makeUniqueEntity<Policy>(createTao())),
      _buffer_matrix(std::move(bufferMatrix)) {
  setType(type);
  Policy::handleError(TaoSetFromOptions(_tao.get()));
}

template <class Policy> void TAOSolver<Policy>::setType(const Type type) {
  TaoType taoType = nullptr;
  switch (type) {
  case Type::lmvm: {
    taoType = TAOLMVM;
    break;
  }
  case Type::nls: {
    taoType = TAONLS;
    break;
  }
  case Type::ntr: {
    taoType = TAONTR;
    break;
  }
  case Type::ntl: {
    taoType = TAONTL;
    break;
  }
  case Type::cg: {
    taoType = TAOCG;
    break;
  }
  case Type::tron: {
    taoType = TAOTRON;
    break;
  }
  case Type::owlqn: {
    taoType = TAOOWLQN;
    break;
  }
  case Type::bmrm: {
    taoType = TAOBMRM;
    break;
  }
  case Type::blmvm: {
    taoType = TAOBLMVM;
    break;
  }
  case Type::bqnls: {
    taoType = TAOBQNLS;
    break;
  }
  case Type::bncg: {
    taoType = TAOBNCG;
    break;
  }
  case Type::bnls: {
    taoType = TAOBNLS;
    break;
  }
  case Type::bntr: {
    taoType = TAOBNTR;
    break;
  }
  case Type::bqnkls: {
    taoType = TAOBQNKLS;
    break;
  }
  case Type::bqnktr: {
    taoType = TAOBQNKTR;
    break;
  }
  case Type::bqnktl: {
    taoType = TAOBQNKTL;
    break;
  }
  case Type::bqpip: {
    taoType = TAOBQPIP;
    break;
  }
  case Type::gpcg: {
    taoType = TAOGPCG;
    break;
  }
  case Type::nm: {
    taoType = TAONM;
    break;
  }
  case Type::pounders: {
    taoType = TAOPOUNDERS;
    break;
  }
  case Type::brgn: {
    taoType = TAOBRGN;
    break;
  }
  case Type::lcl: {
    taoType = TAOLCL;
    break;
  }
  case Type::ssils: {
    taoType = TAOSSILS;
    break;
  }
  case Type::ssfls: {
    taoType = TAOSSFLS;
    break;
  }
  case Type::asils: {
    taoType = TAOASILS;
    break;
  }
  case Type::asfls: {
    taoType = TAOASFLS;
    break;
  }
  case Type::ipm: {
    taoType = TAOIPM;
    break;
  }
  case Type::pdipm: {
    taoType = TAOPDIPM;
    break;
  }
  case Type::shell: {
    taoType = TAOSHELL;
    break;
  }
  case Type::admm: {
    taoType = TAOADMM;
    break;
  }
  case Type::almm: {
    taoType = TAOALMM;
    break;
  }
  }
  Policy::handleError(TaoSetType(_tao.get(), taoType));
}

template <class Policy>
void TAOSolver<Policy>::setBounds(const distributed<vector_type> &lowerBounds,
                                  const distributed<vector_type> &upperBounds) {
  Policy::handleError(TaoSetVariableBounds(
      _tao.get(), lowerBounds.unwrap().data(), upperBounds.unwrap().data()));
}

template <class Policy>
void TAOSolver<Policy>::setConstraints(EqualityConstraints constraints) {
  _constraints.emplace(std::move(constraints));

  Policy::handleError(TaoSetEqualityConstraintsRoutine(
      _tao.get(), _constraints->buffer.first.unwrap().data(),
      &TAOSolver::gradientAdapter, &_constraints->residual));
  Policy::handleError(TaoSetJacobianEqualityRoutine(
      _tao.get(), _constraints->buffer.second.data(),
      _constraints->buffer.second.data(), &TAOSolver::hessianAdapter,
      &_constraints->jacobian));
}

template <class Policy> void TAOSolver<Policy>::checkConvergence() const {
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
distributed<typename TAOSolver<Policy>::vector_type>
TAOSolver<Policy>::solve(ObjectiveFunctor objective, GradientFunctor gradient,
                         HessianFunctor hessian,
                         distributed<vector_type> guess) const {
  Policy::handleError(TaoSetObjectiveRoutine(
      _tao.get(), &TAOSolver::objectiveAdapter, &objective));
  Policy::handleError(TaoSetGradientRoutine(
      _tao.get(), &TAOSolver::gradientAdapter, &gradient));
  Policy::handleError(TaoSetHessianRoutine(
      _tao.get(), _buffer_matrix.data(), _buffer_matrix.data(),
      &TAOSolver::hessianAdapter, &hessian));
  Policy::handleError(TaoSetInitialVector(_tao.get(), guess.unwrap().data()));
  Policy::handleError(TaoSolve(_tao.get()));

  checkConvergence();

  return guess;
}

template <class Policy>
PetscErrorCode TAOSolver<Policy>::objectiveAdapter(Tao, Vec input,
                                                   real_type *output,
                                                   void *context) {
  const auto functionPtr = static_cast<ObjectiveFunctor *>(context);
  (*functionPtr)(distributed<vector_type>(UniqueEntity<Vec>(input, [](Vec) {})),
                 output);
  return 0;
}

template <class Policy>
PetscErrorCode TAOSolver<Policy>::gradientAdapter(Tao, Vec input, Vec output,
                                                  void *context) {
  const auto functionPtr = static_cast<GradientFunctor *>(context);
  distributed<vector_type> wrappedOutput(UniqueEntity<Vec>(output, [](Vec) {}));
  (*functionPtr)(distributed<vector_type>(UniqueEntity<Vec>(input, [](Vec) {})),
                 &wrappedOutput);
  return 0;
}

template <class Policy>
PetscErrorCode TAOSolver<Policy>::hessianAdapter(Tao, Vec input, Mat,
                                                 Mat output, void *context) {
  const auto jacobianPtr = static_cast<HessianFunctor *>(context);
  matrix_type wrappedOutput(UniqueEntity<Mat>(output, [](Mat) {}));
  (*jacobianPtr)(distributed<vector_type>(UniqueEntity<Vec>(input, [](Vec) {})),
                 &wrappedOutput);
  return 0;
}
} // namespace cpppetsc
} // namespace ae108

#endif