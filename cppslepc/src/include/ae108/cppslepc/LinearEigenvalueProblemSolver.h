// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppslepc/EigenPair.h"
#include "ae108/cppslepc/EigenvalueProblemSolverDivergedException.h"
#include <slepceps.h>

namespace ae108 {
namespace cppslepc {

template <class Policy> class LinearEigenvalueProblemSolver {
public:
  using size_type = PetscInt;
  using value_type = PetscScalar;
  using real_type = PetscReal;
  using complex_type = PetscComplex;
  using vector_type = cpppetsc::distributed<cpppetsc::Vector<Policy>>;
  using matrix_type = cpppetsc::Matrix<Policy>;

  explicit LinearEigenvalueProblemSolver();

  enum class Type {
    hermitian = EPS_HEP,
    nonhermitian = EPS_NHEP,
    generalized_hermitian = EPS_GHEP,
    generalized_nonhermitian = EPS_GNHEP,
    generalized_nonhermitian_spd = EPS_PGNHEP,
    generalized_indefinite = EPS_GHIEP
  };

  /**
   * @brief Sets the matrices associated with a standard eigenvalue problem,
   * i.e. Ax = lambda * x.
   * @throw InvalidProblemTypeException if the type does not represent a
   * nongeneralized problem
   */
  void setOperators(const matrix_type *A, const Type type = Type::nonhermitian);

  /**
   * @brief Sets the matrices associated with a generalized eigenvalue
   * problem, i.e. Ax = lambda * Bx.
   * @throw InvalidProblemTypeException if the type does not represent a
   * generalized problem
   */
  void setOperators(const matrix_type *A, const matrix_type *B,
                    const Type type = Type::generalized_nonhermitian);

  /**
   * @brief Solve linear eigenvalue problem.
   * @throw NoOperatorsSetException if the operators have not been set before
   * the call to solve()
   */
  void solve();

  /**
   * @brief Number of discovered eigenpairs.
   */
  size_type numberOfEigenpairs() const;

  /**
   * @brief Get nth eigenpair.
   * @throw InvalidEigenvalueIndexException if `n` is out of range
   */
  void getEigenpair(const size_type n, EigenPair<Policy> *out) const;

  /**
   * @brief Get nth eigenvalue.
   * @throw InvalidEigenvalueIndexException if `n` is out of range
   */
  complex_type getEigenvalue(const size_type n) const;

  /**
   * @brief Returns the internal solver.
   */
  EPS data() const;

private:
  enum class GeneralizedProblem {
    no = 0,
    yes = 1,
  };

  /**
   * @brief Sets the type of the eigenvalue problem.
   * @throw InvalidProblemTypeException if the type does not match the state of
   * `generalized`
   */
  void setType(const Type type, const GeneralizedProblem generalized);

  cpppetsc::UniqueEntity<EPS> eps_;
};

extern template class LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>;
extern template class LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>;

} // namespace cppslepc
} // namespace ae108

#include "ae108/cppslepc/InvalidEigenvalueIndexException.h"
#include "ae108/cppslepc/InvalidProblemTypeException.h"
#include "ae108/cppslepc/NoOperatorsSetException.h"

namespace ae108 {
namespace cppslepc {

template <class Policy>
LinearEigenvalueProblemSolver<Policy>::LinearEigenvalueProblemSolver()
    : eps_([]() {
        auto solver = EPS{};
        Policy::handleError(EPSCreate(Policy::communicator(), &solver));
        return cpppetsc::UniqueEntity<EPS>(
            solver, [](EPS eps) { Policy::handleError(EPSDestroy(&eps)); });
      }()) {
  Policy::handleError(EPSSetFromOptions(this->data()));
}

template <class Policy>
void LinearEigenvalueProblemSolver<Policy>::setOperators(const matrix_type *A,
                                                         const Type type) {
  Policy::handleError(EPSSetOperators(this->data(), A->data(), NULL));
  setType(type, GeneralizedProblem::no);
}

template <class Policy>
void LinearEigenvalueProblemSolver<Policy>::setOperators(const matrix_type *A,
                                                         const matrix_type *B,
                                                         const Type type) {
  Policy::handleError(EPSSetOperators(this->data(), A->data(), B->data()));
  setType(type, GeneralizedProblem::yes);
}

template <class Policy> void LinearEigenvalueProblemSolver<Policy>::solve() {
  const auto numberOfMatrices = [&]() {
    auto st = ST();
    Policy::handleError(EPSGetST(data(), &st));
    auto n = size_type();
    Policy::handleError(STGetNumMatrices(st, &n));
    return n;
  };

  if (numberOfMatrices() == 0) {
    throw NoOperatorsSetException();
  }

  Policy::handleError(EPSSolve(this->data()));

  const auto hasError = [&]() {
    auto reason = EPSConvergedReason{};
    Policy::handleError(EPSGetConvergedReason(this->data(), &reason));
    return reason < 0;
  }();

  if (hasError) {
    throw EigenvalueProblemSolverDivergedException{};
  }
}

template <class Policy>
typename LinearEigenvalueProblemSolver<Policy>::size_type
LinearEigenvalueProblemSolver<Policy>::numberOfEigenpairs() const {
  auto size = size_type{};
  Policy::handleError(EPSGetConverged(this->data(), &size));
  return size;
}

template <class Policy>
void LinearEigenvalueProblemSolver<Policy>::setType(
    const Type type, const GeneralizedProblem generalized) {
  const auto generalizedType =
      bool{type == Type::generalized_hermitian ||
           type == Type::generalized_nonhermitian ||
           type == Type::generalized_nonhermitian_spd ||
           type == Type::generalized_indefinite};
  if (generalizedType != (generalized == GeneralizedProblem::yes)) {
    throw InvalidProblemTypeException();
  }
  Policy::handleError(
      EPSSetProblemType(this->data(), static_cast<EPSProblemType>(type)));
}

template <class Policy>
void LinearEigenvalueProblemSolver<Policy>::getEigenpair(
    const size_type n, EigenPair<Policy> *out) const {
  if (n < 0 || n >= numberOfEigenpairs()) {
    throw InvalidEigenvalueIndexException();
  }

#ifdef AE108_PETSC_COMPLEX
  Policy::handleError(EPSGetEigenpair(this->data(), n, &out->value, NULL,
                                      out->vector.unwrap().data(), NULL));
#else
  value_type real, imag;
  Policy::handleError(EPSGetEigenpair(this->data(), n, &real, &imag,
                                      out->vector_real.unwrap().data(),
                                      out->vector_imag.unwrap().data()));
  out->value.real(real);
  out->value.imag(imag);
#endif
}

template <class Policy>
typename LinearEigenvalueProblemSolver<Policy>::complex_type
LinearEigenvalueProblemSolver<Policy>::getEigenvalue(const size_type n) const {
  if (n < 0 || n >= numberOfEigenpairs()) {
    throw InvalidEigenvalueIndexException();
  }

  auto eigenvalue = complex_type();

#ifdef AE108_PETSC_COMPLEX
  Policy::handleError(EPSGetEigenvalue(this->data(), n, &eigenvalue, NULL));
#else
  value_type real, imag;
  Policy::handleError(EPSGetEigenvalue(this->data(), n, &real, &imag));
  eigenvalue.real(real);
  eigenvalue.imag(imag);
#endif

  return eigenvalue;
}

template <class Policy>
EPS LinearEigenvalueProblemSolver<Policy>::data() const {
  return eps_.get();
}

} // namespace cppslepc
} // namespace ae108
