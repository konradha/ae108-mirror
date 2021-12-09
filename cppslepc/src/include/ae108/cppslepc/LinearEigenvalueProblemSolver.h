// © 2021 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppslepc/EigenPair.h"
#include "ae108/cppslepc/EigenvalueProblemSolverDivergedException.h"
#include <slepc/slepceps.h>

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

  /**
   * @brief Sets the matrices associated with a standard eigenvalue problem,
   * i.e. Ax = lambda * x.
   */
  void setOperators(const matrix_type *A);

  /**
   * @brief Sets the matrices associated with a generalized eigenvalue
   * problem, i.e. Ax = lambda * Bx.
   */
  void setOperators(const matrix_type *A, const matrix_type *B);

  /**
   * @brief Solve linear eigenvalue problem and returns the number of converged
   * eigenpairs.
   */
  size_type solve() const;

  /**
   * @brief Get nth eigenpair.
   */
  void getEigenpair(const size_type n, EigenPair<Policy> *out) const;

  /**
   * @brief Get nth eigenvalue.
   */
  complex_type getEigenvalue(const size_type n) const;

  /**
   * @brief Returns the internal solver.
   */
  EPS data() const;

private:
  cpppetsc::UniqueEntity<EPS> eps_;
};

extern template class LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>;
extern template class LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>;

} // namespace cppslepc
} // namespace ae108

namespace ae108 {
namespace cppslepc {

template <class Policy>
LinearEigenvalueProblemSolver<Policy>::LinearEigenvalueProblemSolver()
    : eps_([]() {
        auto solver = EPS{};
        Policy::handleError(EPSCreate(Policy::communicator(), &solver));
        return cpppetsc::UniqueEntity<EPS>(
            solver, [](EPS eps) { Policy::handleError(EPSDestroy(&eps)); });
      }()) {}

template <class Policy>
void LinearEigenvalueProblemSolver<Policy>::setOperators(
    const typename LinearEigenvalueProblemSolver<Policy>::matrix_type *A) {
  Policy::handleError(EPSSetOperators(this->data(), A->data(), NULL));
}

template <class Policy>
void LinearEigenvalueProblemSolver<Policy>::setOperators(
    const typename LinearEigenvalueProblemSolver<Policy>::matrix_type *A,
    const typename LinearEigenvalueProblemSolver<Policy>::matrix_type *B) {
  Policy::handleError(EPSSetOperators(this->data(), A->data(), B->data()));
}

template <class Policy>
typename LinearEigenvalueProblemSolver<Policy>::size_type
LinearEigenvalueProblemSolver<Policy>::solve() const {

  Policy::handleError(EPSSetFromOptions(this->data()));

  Policy::handleError(EPSSolve(this->data()));

  const auto hasError = [&]() {
    auto reason = EPSConvergedReason{};
    Policy::handleError(EPSGetConvergedReason(this->data(), &reason));
    return reason < 0;
  }();

  if (hasError) {
    throw EigenvalueProblemSolverDivergedException{};
  }

  return [&]() {
    auto size = size_type{};
    Policy::handleError(EPSGetConverged(this->data(), &size));
    return size;
  }();
}

template <class Policy>
void LinearEigenvalueProblemSolver<Policy>::getEigenpair(
    const size_type n, EigenPair<Policy> *out) const {

#ifdef AE108_PETSC_COMPLEX
  EPSGetEigenpair(this->data(), n, &out->value, NULL,
                  out->vector.unwrap().data(), NULL);
#else
  value_type real, imag;
  EPSGetEigenpair(this->data(), n, &real, &imag,
                  out->vector_real.unwrap().data(),
                  out->vector_imag.unwrap().data());
  out->value.real(real);
  out->value.imag(imag);
#endif
}

template <class Policy>
typename LinearEigenvalueProblemSolver<Policy>::complex_type
LinearEigenvalueProblemSolver<Policy>::getEigenvalue(const size_type n) const {

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