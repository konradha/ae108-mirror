#pragma once

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include <complex>
#include <slepc/slepceps.h>
#include <vector>

namespace ae108 {
namespace cppslepc {

template <class Policy> class LinearEigenvalueProblemSolver {
public:
  using size_type = PetscInt;
  using value_type = PetscScalar;
  using real_type = PetscReal;

  explicit LinearEigenvalueProblemSolver();

  /**
   * @brief Returns the detected eigenvalues of the generalized eigenvalue
   * problem A x = lambda * B x.
   */
  std::vector<std::complex<real_type>>
  solve(const cpppetsc::Matrix<Policy> &A,
        const cpppetsc::Matrix<Policy> &B) const;

  /**
   * @brief Returns the detected eigenvalues of A.
   */
  std::vector<std::complex<real_type>>
  solve(const cpppetsc::Matrix<Policy> &A) const;

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

#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

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
  Policy::handleError(EPSSetFromOptions(eps_.get()));
}

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
LinearEigenvalueProblemSolver<Policy>::solve(
    const cpppetsc::Matrix<Policy> &A,
    const cpppetsc::Matrix<Policy> &B) const {
  Policy::handleError(EPSSetOperators(eps_.get(), A.data(), B.data()));
  Policy::handleError(EPSSolve(eps_.get()));

  const auto size = [&]() {
    auto size = size_type{};
    Policy::handleError(EPSGetConverged(eps_.get(), &size));
    return size;
  }();

  namespace rv = ranges::cpp20::views;
  return rv::iota(0, size) | rv::transform([&](const size_type index) {
           auto eigenvalue = std::pair<PetscScalar, PetscScalar>();
           Policy::handleError(EPSGetEigenvalue(
               eps_.get(), index, &eigenvalue.first, &eigenvalue.second));
           return
#ifdef AE108_PETSC_COMPLEX
               eigenvalue.first
#else
               std::complex<real_type> {
             eigenvalue.first, eigenvalue.second
           }
#endif
               ;
         }) |
         ranges::to<std::vector>();
}

template <class Policy>
std::vector<
    std::complex<typename LinearEigenvalueProblemSolver<Policy>::real_type>>
LinearEigenvalueProblemSolver<Policy>::solve(
    const cpppetsc::Matrix<Policy> &A) const {
  Policy::handleError(EPSSetOperators(eps_.get(), A.data(), nullptr));
  Policy::handleError(EPSSolve(eps_.get()));

  const auto size = [&]() {
    auto size = size_type{};
    Policy::handleError(EPSGetConverged(eps_.get(), &size));
    return size;
  }();

  namespace rv = ranges::cpp20::views;
  return rv::iota(0, size) | rv::transform([&](const size_type index) {
           auto eigenvalue = std::pair<PetscScalar, PetscScalar>();
           Policy::handleError(EPSGetEigenvalue(
               eps_.get(), index, &eigenvalue.first, &eigenvalue.second));
           return
#ifdef AE108_PETSC_COMPLEX
               eigenvalue.first
#else
               std::complex<real_type> {
             eigenvalue.first, eigenvalue.second
           }
#endif
               ;
         }) |
         ranges::to<std::vector>();
}

template <class Policy>
EPS LinearEigenvalueProblemSolver<Policy>::data() const {
  return eps_.get();
}

} // namespace cppslepc
} // namespace ae108