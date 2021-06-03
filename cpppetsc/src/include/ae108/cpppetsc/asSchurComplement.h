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
 * @param indices Indices.
 */
template <class Policy>
Matrix<Policy> asSchurComplement(
    const Matrix<Policy> *matrix,
    const std::vector<typename Matrix<Policy>::size_type> indices);

extern template Matrix<SequentialComputePolicy> asSchurComplement(
    const Matrix<SequentialComputePolicy> *matrix,
    const std::vector<typename Matrix<SequentialComputePolicy>::size_type>
        indices);

extern template Matrix<ParallelComputePolicy> asSchurComplement(
    const Matrix<ParallelComputePolicy> *matrix,
    const std::vector<typename Matrix<ParallelComputePolicy>::size_type>
        indices);

} // namespace cpppetsc
} // namespace ae108

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
Matrix<Policy>
asSchurComplement(const Matrix<Policy> *matrix,
                  std::vector<typename Matrix<Policy>::size_type> indices) {
  assert(matrix);

  auto is0 = IS();
  Policy::handleError(ISCreateGeneral(PETSC_COMM_WORLD, indices.size(),
                                      indices.data(), PETSC_COPY_VALUES, &is0));

  auto is1 = IS();
  Policy::handleError(ISComplement(is0, 0, matrix->size().first, &is1));

  auto mat = Mat{};
  Policy::handleError(MatGetSchurComplement(
      matrix->data(), is0, is0, is1, is1, MAT_INITIAL_MATRIX, &mat,
      MAT_SCHUR_COMPLEMENT_AINV_DIAG, MAT_IGNORE_MATRIX, NULL));

  return Matrix<Policy>(makeUniqueEntity<Policy>(mat));
}

} // namespace cpppetsc
} // namespace ae108