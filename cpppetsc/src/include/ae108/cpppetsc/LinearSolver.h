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

#include "ae108/cpppetsc/LinearSolverDivergedException.h"
#include "ae108/cpppetsc/Matrix_fwd.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector_fwd.h"
#include <memory>
#include <petscksp.h>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief This class can be used to solve linear systems Ax = b.
 */
template <class Policy> class LinearSolver {
public:
  using matrix_type = Matrix<Policy>;
  using vector_type = Vector<Policy>;

  /**
   * @brief Initialize the solver with a matrix A.
   *
   * @remark The solver keeps a reference to that matrix.
   * @param mat A pointer to the matrix A.
   */
  explicit LinearSolver(const matrix_type *mat);

  /**
   * @brief Calls solve with a suitable result vector. Passes on
   * the rest of the arguments.
   *
   * @param rhs The vector b.
   * @return The solution of the linear system.
   * @throw LinearSolverDivergedException if the solve diverged
   */
  distributed<vector_type> solve(const distributed<vector_type> &rhs) const;

  /**
   * @brief Solve the linear system Ax = b.
   *
   * @param rhs The vector b.
   * @param result Used to store the result (instead of creating a new vector).
   * @return The solution of the linear system.
   * @throw LinearSolverDivergedException if the solve diverged
   */
  distributed<vector_type> solve(const distributed<vector_type> &rhs,
                                 distributed<vector_type> result) const;

private:
  static KSP createKSP();

  UniqueEntity<KSP> _ksp;
};

extern template class LinearSolver<SequentialComputePolicy>;
extern template class LinearSolver<ParallelComputePolicy>;
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

template <class Policy> KSP LinearSolver<Policy>::createKSP() {
  KSP ksp;
  Policy::handleError(KSPCreate(Policy::communicator(), &ksp));
  return ksp;
}

template <class Policy>
LinearSolver<Policy>::LinearSolver(const matrix_type *mat)
    : _ksp(makeUniqueEntity<Policy>(createKSP())) {
  Policy::handleError(KSPSetFromOptions(_ksp.get()));

  // the sizes of all matrices passed to KSP must be the same, so we only permit
  // setting them once in the constructor
  Policy::handleError(KSPSetOperators(_ksp.get(), mat->data(), mat->data()));
}

template <class Policy>
distributed<typename LinearSolver<Policy>::vector_type>
LinearSolver<Policy>::solve(const distributed<vector_type> &rhs) const {
  return solve(rhs,
               tag<DistributedTag>(vector_type::fromLayoutOf(rhs.unwrap())));
}

template <class Policy>
distributed<typename LinearSolver<Policy>::vector_type>
LinearSolver<Policy>::solve(const distributed<vector_type> &rhs,
                            distributed<vector_type> result) const {
  Policy::handleError(
      KSPSolve(_ksp.get(), rhs.unwrap().data(), result.unwrap().data()));

  const auto errorCode = [this]() {
    auto code = KSPConvergedReason{};
    Policy::handleError(KSPGetConvergedReason(_ksp.get(), &code));
    return code;
  }();

  if (errorCode < 0) {
    throw LinearSolverDivergedException{};
  }
  return result;
}
} // namespace cpppetsc
} // namespace ae108