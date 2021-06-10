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
#include <algorithm>
#include <cassert>
#include <petscksp.h>

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

  auto wrappedIn =
      distributed<Vector<Policy>>(UniqueEntity<Vec>(in, [](Vec) {}));
  auto wrappedOut =
      distributed<Vector<Policy>>(UniqueEntity<Vec>(out, [](Vec) {}));

  data->solver.solve(
      distributed<Vector<Policy>>(UniqueEntity<Vec>(in, [](Vec) {})),
      distributed<Vector<Policy>>(UniqueEntity<Vec>(out, [](Vec) {}))
  );

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