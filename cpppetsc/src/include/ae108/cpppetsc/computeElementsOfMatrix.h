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
 * @brief Computes the elements of the given matrix by evaluating matrix-vector
 * products.
 */
template <class Policy>
Matrix<Policy> computeElementsOfMatrix(const Matrix<Policy> &matrix);

extern template Matrix<SequentialComputePolicy>
computeElementsOfMatrix(const Matrix<SequentialComputePolicy> &);

extern template Matrix<ParallelComputePolicy>
computeElementsOfMatrix(const Matrix<ParallelComputePolicy> &);

} // namespace cpppetsc
} // namespace ae108

namespace ae108 {
namespace cpppetsc {

template <class Policy>
cpppetsc::Matrix<Policy>
computeElementsOfMatrix(const cpppetsc::Matrix<Policy> &matrix) {
  auto result = Mat{};
  Policy::handleError(MatComputeOperator(matrix.data(), MATAIJ, &result));
  return cpppetsc::Matrix<Policy>(makeUniqueEntity<Policy>(result));
}

} // namespace cpppetsc
} // namespace ae108