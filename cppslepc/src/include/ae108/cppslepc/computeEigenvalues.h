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
computeEigenvalues(const cpppetsc::Matrix<Policy> &A,
                   const cpppetsc::Matrix<Policy> &B);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
                   const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &);

extern template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &,
                   const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &);

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
computeEigenvalues(const cpppetsc::Matrix<Policy> &A) {
  auto solver = LinearEigenvalueProblemSolver<Policy>{};

  solver.setOperators(&A);

  return detail::solve(std::move(solver));
}

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<Policy> &A,
                   const cpppetsc::Matrix<Policy> &B) {
  auto solver = LinearEigenvalueProblemSolver<Policy>{};

  solver.setOperators(&A, &B);

  return detail::solve(std::move(solver));
}

} // namespace cppslepc
} // namespace ae108