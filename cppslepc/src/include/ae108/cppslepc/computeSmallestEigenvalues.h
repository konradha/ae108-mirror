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
#include "ae108/cppslepc/LinearEigenvalueProblemSolver.h"
#include <complex>
#include <vector>

namespace ae108 {
namespace cppslepc {

/**
 * @brief Returns the n smallest found eigenvalues (by magnitude) of the
 * eigenvalue problem A x = lambda * x.
 * @remark The selection of returned eigenvalues may depend on the
 * `numberOfEigenvalues` parameter.
 *
 * @throw EigenvalueProblemSolverDivergedException if the solver failed to find
 * a solution.
 */
template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeSmallestEigenvalues(const cpppetsc::Matrix<Policy> &A,
                           const std::size_t numberOfEigenvalues);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>::real_type>>
computeSmallestEigenvalues(
    const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
    const std::size_t numberOfEigenvalues);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeSmallestEigenvalues(
    const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &,
    const std::size_t numberOfEigenvalues);

/**
 * @brief Returns the n smallest found eigenvalues (by magnitude) of the
 * generalized eigenvalue problem A x = lambda * B x.
 * @remark The selection of returned eigenvalues may depend on the
 * `numberOfEigenvalues` parameter.
 *
 * @throw EigenvalueProblemSolverDivergedException if the solver failed to find
 * a solution.
 */
template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeSmallestEigenvalues(const cpppetsc::Matrix<Policy> &A,
                           const cpppetsc::Matrix<Policy> &B,
                           const std::size_t numberOfEigenvalues);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>::real_type>>
computeSmallestEigenvalues(
    const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
    const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
    const std::size_t numberOfEigenvalues);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeSmallestEigenvalues(
    const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &,
    const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &,
    const std::size_t numberOfEigenvalues);

} // namespace cppslepc
} // namespace ae108

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <slepceps.h>

namespace ae108 {
namespace cppslepc {

namespace detail {

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
solve(LinearEigenvalueProblemSolver<Policy> solver) {
  using solver_type = LinearEigenvalueProblemSolver<Policy>;
  using size_type = typename solver_type::size_type;

  namespace rv = ranges::cpp20::views;
  solver.solve();
  return rv::iota(0, solver.numberOfEigenpairs()) |
         rv::transform([&](const size_type index) {
           return solver.getEigenvalue(index);
         }) |
         ranges::to<std::vector>();
}
} // namespace detail

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeSmallestEigenvalues(const cpppetsc::Matrix<Policy> &A,
                           const std::size_t numberOfEigenvalues) {
  if (!numberOfEigenvalues)
    return {};

  auto solver = LinearEigenvalueProblemSolver<Policy>{};
  Policy::handleError(
      EPSSetWhichEigenpairs(solver.data(), EPS_SMALLEST_MAGNITUDE));
  Policy::handleError(EPSSetDimensions(solver.data(), numberOfEigenvalues,
                                       PETSC_DEFAULT, PETSC_DEFAULT));

  solver.setOperators(&A);

  return detail::solve(std::move(solver));
}

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeSmallestEigenvalues(const cpppetsc::Matrix<Policy> &A,
                           const cpppetsc::Matrix<Policy> &B,
                           const std::size_t numberOfEigenvalues) {
  if (!numberOfEigenvalues)
    return {};

  auto solver = LinearEigenvalueProblemSolver<Policy>{};
  Policy::handleError(
      EPSSetWhichEigenpairs(solver.data(), EPS_SMALLEST_MAGNITUDE));
  Policy::handleError(EPSSetDimensions(solver.data(), numberOfEigenvalues,
                                       PETSC_DEFAULT, PETSC_DEFAULT));

  solver.setOperators(&A, &B);

  return detail::solve(std::move(solver));
}

} // namespace cppslepc
} // namespace ae108