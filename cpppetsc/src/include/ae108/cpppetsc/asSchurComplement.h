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
 * @brief Returns a matrix that behaves like the Schur complement:
 * M_11 - M_10 * M_00^-1 * M_01
 *
 * Note that the matrix is not actually computed, but the returned matrix
 * stores a reference to the given matrices to compute the result of
 * operations on demand.
 *
 * @param matrix_00 Valid nonzero pointer.
 * @param matrix_01 Valid nonzero pointer.
 * @param matrix_10 Valid nonzero pointer.
 * @param matrix_11 Valid nonzero pointer.
 */
template <class Policy>
Matrix<Policy> asSchurComplement(const Matrix<Policy> *matrix_00,
                                 const Matrix<Policy> *matrix_01,
                                 const Matrix<Policy> *matrix_10,
                                 const Matrix<Policy> *matrix_11);

extern template Matrix<SequentialComputePolicy>
asSchurComplement(const Matrix<SequentialComputePolicy> *matrix_00,
                  const Matrix<SequentialComputePolicy> *matrix_01,
                  const Matrix<SequentialComputePolicy> *matrix_10,
                  const Matrix<SequentialComputePolicy> *matrix_11);

extern template Matrix<ParallelComputePolicy>
asSchurComplement(const Matrix<ParallelComputePolicy> *matrix_00,
                  const Matrix<ParallelComputePolicy> *matrix_01,
                  const Matrix<ParallelComputePolicy> *matrix_10,
                  const Matrix<ParallelComputePolicy> *matrix_11);

/**
 * @brief Returns a matrix that behaves like the Schur complement:
 * M_11 - M_10 * M_00^-1 * M_01
 *
 * Note that the matrix is not actually computed, but the returned matrix
 * stores a reference to the given matrices to compute the result of
 * operations on demand.
 *
 * @param matrix Valid nonzero pointer.
 * @param indices The row/column indices that define the matrix M_11.
 * Only indices in the local row range of the matrix need to be provided.
 */
template <class Policy>
Matrix<Policy> asSchurComplement(
    const Matrix<Policy> *matrix,
    const std::vector<typename Matrix<Policy>::size_type> &indices);

extern template Matrix<SequentialComputePolicy> asSchurComplement(
    const Matrix<SequentialComputePolicy> *,
    const std::vector<typename Matrix<SequentialComputePolicy>::size_type> &);

extern template Matrix<ParallelComputePolicy> asSchurComplement(
    const Matrix<ParallelComputePolicy> *,
    const std::vector<typename Matrix<ParallelComputePolicy>::size_type> &);

} // namespace cpppetsc
} // namespace ae108

#include <algorithm>
#include <cassert>
#include <petscksp.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy>
Matrix<Policy> asSchurComplement(const Matrix<Policy> *matrix_00,
                                 const Matrix<Policy> *matrix_01,
                                 const Matrix<Policy> *matrix_10,
                                 const Matrix<Policy> *matrix_11) {
  assert(matrix_00);
  assert(matrix_01);
  assert(matrix_10);
  assert(matrix_11);

  auto mat = Mat{};
  Policy::handleError(MatCreateSchurComplement(
      matrix_00->data(), matrix_00->data(), matrix_01->data(),
      matrix_10->data(), matrix_11->data(), &mat));
  return Matrix<Policy>(makeUniqueEntity<Policy>(mat));
}

template <class Policy>
Matrix<Policy> asSchurComplement(
    const Matrix<Policy> *matrix,
    const std::vector<typename Matrix<Policy>::size_type> &indices) {
  assert(matrix);

  const auto localRowRange = matrix->localRowRange();

  const auto indices_11 = [&]() {
    using size_type = typename Matrix<Policy>::size_type;

    const auto isLocal = [&](const typename Matrix<Policy>::size_type index) {
      return localRowRange.first <= index && index < localRowRange.second;
    };
    const auto localSize =
        std::count_if(indices.begin(), indices.end(), isLocal);

    auto is = IS{};
    auto localIndices = static_cast<size_type *>(nullptr);
    Policy::handleError(
        PetscMalloc(sizeof(size_type) * localSize, &localIndices));
    std::copy_if(indices.begin(), indices.end(), localIndices, isLocal);
    Policy::handleError(ISCreateGeneral(Policy::communicator(), localSize,
                                        localIndices, PETSC_OWN_POINTER, &is));

    return makeUniqueEntity<Policy>(is);
  }();

  const auto indices_00 = [&]() {
    auto is = IS{};
    Policy::handleError(ISComplement(indices_11.get(), localRowRange.first,
                                     localRowRange.second, &is));
    return makeUniqueEntity<Policy>(is);
  }();

  auto mat = Mat{};
  Policy::handleError(MatGetSchurComplement(
      matrix->data(), indices_00.get(), indices_00.get(), indices_11.get(),
      indices_11.get(), MAT_INITIAL_MATRIX, &mat,
      MAT_SCHUR_COMPLEMENT_AINV_DIAG, MAT_IGNORE_MATRIX, nullptr));

  return Matrix<Policy>(makeUniqueEntity<Policy>(mat));
}

} // namespace cpppetsc
} // namespace ae108