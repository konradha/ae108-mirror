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

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Returns a matrix that behaves like the inverse of the provided matrix.
 *
 * Note that the inverse matrix is not actually computed, but the returned
 * matrix stores a reference to the given matrix to compute the result of
 * operations on demand.
 *
 * @param matrix Valid nonzero pointer.
 */
template <class Policy>
Matrix<Policy> asInverseMatrix(const Matrix<Policy> *matrix);

extern template Matrix<SequentialComputePolicy>
asInverseMatrix(const Matrix<SequentialComputePolicy> *);

extern template Matrix<ParallelComputePolicy>
asInverseMatrix(const Matrix<ParallelComputePolicy> *);

} // namespace cpppetsc
} // namespace ae108

#include "ae108/cpppetsc/LinearSolver.h"
#include <cassert>
#include <memory>
#include <petscksp.h>
#include <utility>

namespace ae108 {
namespace cpppetsc {

namespace detail {
namespace asInverseMatrix {

/**
 * @brief Stores the data needed to compute the matrix by vector
 * product.
 */
template <class Policy> struct Data { LinearSolver<Policy> solver; };

/**
 * @brief Computes `matrix^-1 * in` using the solver in the context.
 */
template <class Policy> PetscErrorCode multiply(Mat mat, Vec in, Vec out) {
  const auto data = [&]() {
    auto context = static_cast<Data<Policy> *>(nullptr);
    Policy::handleError(MatShellGetContext(mat, &context));
    return context;
  }();

  data->solver.solve(
      distributed<Vector<Policy>>(UniqueEntity<Vec>(in, [](Vec) {})),
      distributed<Vector<Policy>>(UniqueEntity<Vec>(out, [](Vec) {})));

  return PetscErrorCode{};
}

/**
 * @brief Deletes the context.
 */
template <class Policy> PetscErrorCode destroy(Mat mat) {
  const auto data = [&]() {
    auto context = static_cast<Data<Policy> *>(nullptr);
    Policy::handleError(MatShellGetContext(mat, &context));
    return context;
  }();

  delete data;

  return PetscErrorCode{};
}

} // namespace asInverseMatrix
} // namespace detail

template <class Policy>
Matrix<Policy> asInverseMatrix(const Matrix<Policy> *matrix) {
  assert(matrix);

  using size_type = typename Matrix<Policy>::size_type;

  const auto localSize = [&]() {
    auto size = std::pair<size_type, size_type>{};
    Policy::handleError(
        MatGetLocalSize(matrix->data(), &size.first, &size.second));
    return size;
  }();

  auto shell = [&]() {
    auto mat = Mat{};
    auto data = std::unique_ptr<detail::asInverseMatrix::Data<Policy>>(
        new detail::asInverseMatrix::Data<Policy>{
            LinearSolver<Policy>(matrix)});

    Policy::handleError(MatCreateShell(Policy::communicator(), localSize.first,
                                       localSize.second, PETSC_DETERMINE,
                                       PETSC_DETERMINE, data.get(), &mat));

    auto result = Matrix<Policy>(makeUniqueEntity<Policy>(mat));

    Policy::handleError(MatShellSetOperation(
        result.data(), MATOP_DESTROY,
        (void (*)(void))detail::asInverseMatrix::destroy<Policy>));
    data.release();

    return result;
  }();

  Policy::handleError(MatShellSetOperation(
      shell.data(), MATOP_MULT,
      (void (*)(void))detail::asInverseMatrix::multiply<Policy>));

  return shell;
}

} // namespace cpppetsc
} // namespace ae108