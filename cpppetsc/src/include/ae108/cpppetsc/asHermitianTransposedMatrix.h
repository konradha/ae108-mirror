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
 * @brief Returns a matrix that behaves like the Hermitian transpose.
 *
 * Note that the matrix is not actually computed, but the returned matrix
 * stores a reference to the given matrix to compute the result of
 * operations on demand.
 *
 * @param matrix Valid nonzero pointer.
 */
template <class Policy>
Matrix<Policy> asHermitianTransposedMatrix(const Matrix<Policy> *matrix);

extern template Matrix<SequentialComputePolicy>
asHermitianTransposedMatrix(const Matrix<SequentialComputePolicy> *matrix);

extern template Matrix<ParallelComputePolicy>
asHermitianTransposedMatrix(const Matrix<ParallelComputePolicy> *matrix);

} // namespace cpppetsc
} // namespace ae108

#include <algorithm>
#include <cassert>

namespace ae108 {
namespace cpppetsc {

template <class Policy>
Matrix<Policy> asHermitianTransposedMatrix(const Matrix<Policy> *matrix) {
  assert(matrix);

  auto mat = Mat{};
  Policy::handleError(MatCreateHermitianTranspose(matrix->data(), &mat));
  return Matrix<Policy>(makeUniqueEntity<Policy>(mat));
}

} // namespace cpppetsc
} // namespace ae108