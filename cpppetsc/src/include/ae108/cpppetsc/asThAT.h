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
 * @brief Returns a matrix that behaves like the transformed matrix
 * T^H * A * T, where ^H denotes the Hermitian transpose.
 *
 * Note that the matrix is not actually computed, but the return matrix
 * stores a reference to the given matrices to compute the result of
 * operations on demand.
 *
 * @param matrix Valid nonzero pointer.
 * @param transform Valid nonzero pointer.
 */
template <class Policy>
Matrix<Policy> asThAT(const Matrix<Policy> *matrix,
                      const Matrix<Policy> *transform);

extern template Matrix<SequentialComputePolicy>
asThAT(const Matrix<SequentialComputePolicy> *,
       const Matrix<SequentialComputePolicy> *);
extern template Matrix<ParallelComputePolicy>
asThAT(const Matrix<ParallelComputePolicy> *,
       const Matrix<ParallelComputePolicy> *);

} // namespace cpppetsc
} // namespace ae108

#include "ae108/cpppetsc/createTransformInput.h"
#include "ae108/cpppetsc/createTransformOutput.h"
#include <cassert>
#include <petscmat.h>
#include <utility>

namespace ae108 {
namespace cpppetsc {

namespace detail {
namespace asThAT {

/**
 * @brief Stores the data needed to compute the matrix by vector
 * product. Also contains two vectors to store intermediate
 * results.
 */
template <class Policy> struct Data {
  Matrix<Policy> matrix;
  Matrix<Policy> transform;

  distributed<Vector<Policy>> x;
  distributed<Vector<Policy>> y;
};

/**
 * @brief Computes `transform^H * matrix * transform * in` using the
 * matrices in the context.
 */
template <class Policy> PetscErrorCode multiply(Mat mat, Vec in, Vec out) {
  const auto data = [&]() {
    auto context = static_cast<Data<Policy> *>(nullptr);
    Policy::handleError(MatShellGetContext(mat, &context));
    return context;
  }();

  Policy::handleError(
      MatMult(data->transform.data(), in, data->x.unwrap().data()));
  Policy::handleError(MatMult(data->matrix.data(), data->x.unwrap().data(),
                              data->y.unwrap().data()));
  Policy::handleError(MatMultHermitianTranspose(data->transform.data(),
                                                data->y.unwrap().data(), out));

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

/**
 * @brief Increases the reference count of the contained `Mat` of the given
 * matrix and returns a matrix wrapping the result.
 */
template <class Policy>
Matrix<Policy> shallowCopy(const Matrix<Policy> *matrix) {
  assert(matrix);
  Policy::handleError(PetscObjectReference((PetscObject)matrix->data()));
  return Matrix<Policy>(makeUniqueEntity<Policy>(matrix->data()));
}

} // namespace asThAT
} // namespace detail

template <class Policy>
Matrix<Policy> asThAT(const Matrix<Policy> *matrix,
                      const Matrix<Policy> *transform) {
  assert(matrix);
  assert(transform);

  using size_type = typename Matrix<Policy>::size_type;

  const auto localSize = [&]() {
    auto size = std::pair<size_type, size_type>{};
    Policy::handleError(
        MatGetLocalSize(transform->data(), &size.first, &size.second));
    return size;
  }();

  auto shell = [&]() {
    auto mat = Mat{};
    auto data = std::unique_ptr<detail::asThAT::Data<Policy>>(
        new detail::asThAT::Data<Policy>{
            detail::asThAT::shallowCopy<Policy>(matrix),
            detail::asThAT::shallowCopy<Policy>(transform),
            createTransformInput(*matrix), createTransformOutput(*matrix)});

    Policy::handleError(MatCreateShell(Policy::communicator(), localSize.second,
                                       localSize.second, PETSC_DETERMINE,
                                       PETSC_DETERMINE, data.get(), &mat));

    auto result = Matrix<Policy>(makeUniqueEntity<Policy>(mat));

    Policy::handleError(
        MatShellSetOperation(result.data(), MATOP_DESTROY,
                             (void (*)(void))detail::asThAT::destroy<Policy>));
    data.release();

    return result;
  }();

  Policy::handleError(
      MatShellSetOperation(shell.data(), MATOP_MULT,
                           (void (*)(void))detail::asThAT::multiply<Policy>));

  return shell;
}

} // namespace cpppetsc
} // namespace ae108