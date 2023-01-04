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