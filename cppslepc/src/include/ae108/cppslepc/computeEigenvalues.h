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
#include "ae108/cppslepc/EigenvalueProblemSolverDivergedException.h"
#include "ae108/cppslepc/LinearEigenvalueProblemSolver.h"
#include <complex>
#include <vector>

namespace ae108 {
namespace cppslepc {

/**
 * @brief Returns the detected eigenvalues of A.
 *
 * @throw EigenvalueProblemSolverDivergedException if the solver failed to find
 * a solution.
 */
template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<Policy> &A);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &);

/**
 * @brief Returns the detected eigenvalues of the generalized eigenvalue
 * problem A x = lambda * B x.
 *
 * @throw EigenvalueProblemSolverDivergedException if the solver failed to find
 * a solution.
 */
template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeGeneralizedEigenvalues(const cpppetsc::Matrix<Policy> &A,
                              const cpppetsc::Matrix<Policy> &B,
                              const std::size_t number_of_eigenvalues);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>::real_type>>
computeGeneralizedEigenvalues(
    const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
    const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
    const std::size_t number_of_eigenvalues);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeGeneralizedEigenvalues(
    const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &,
    const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &,
    const std::size_t number_of_eigenvalues);

} // namespace cppslepc
} // namespace ae108

#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <slepceps.h>

namespace ae108 {
namespace cppslepc {

namespace detail {

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
solve(const LinearEigenvalueProblemSolver<Policy> &solver) {
  using solver_type = LinearEigenvalueProblemSolver<Policy>;
  using size_type = typename solver_type::size_type;

  Policy::handleError(EPSSolve(solver.data()));

  const auto hasError = [&]() {
    auto reason = EPSConvergedReason{};
    Policy::handleError(EPSGetConvergedReason(solver.data(), &reason));
    return reason < 0;
  }();

  if (hasError) {
    throw EigenvalueProblemSolverDivergedException{};
  }

  const auto size = [&]() {
    auto size = size_type{};
    Policy::handleError(EPSGetConverged(solver.data(), &size));
    return size;
  }();

  namespace rv = ranges::cpp20::views;
  return rv::iota(0, size) | rv::transform([&](const size_type index) {
           auto eigenvalue = std::pair<PetscScalar, PetscScalar>();
           Policy::handleError(EPSGetEigenvalue(
               solver.data(), index, &eigenvalue.first, &eigenvalue.second));
           return
#ifdef AE108_PETSC_COMPLEX
               eigenvalue.first
#else
               std::complex<typename solver_type::real_type> {
             eigenvalue.first, eigenvalue.second
           }
#endif
               ;
         }) |
         ranges::to<std::vector>();
}
} // namespace detail

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<Policy> &A) {
  auto solver = LinearEigenvalueProblemSolver<Policy>{};

  Policy::handleError(EPSSetOperators(solver.data(), A.data(), nullptr));
  return detail::solve(solver);
}

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeGeneralizedEigenvalues(const cpppetsc::Matrix<Policy> &A,
                              const cpppetsc::Matrix<Policy> &B,
                              const std::size_t number_of_eigenvalues) {
  auto solver = LinearEigenvalueProblemSolver<Policy>{};

  Policy::handleError(
      EPSSetWhichEigenpairs(solver.data(), EPS_SMALLEST_MAGNITUDE));
  Policy::handleError(EPSSetDimensions(solver.data(), number_of_eigenvalues,
                                       PETSC_DECIDE, PETSC_DECIDE));
  Policy::handleError(EPSSetOperators(solver.data(), A.data(), B.data()));
  return detail::solve(solver);
}

} // namespace cppslepc
} // namespace ae108