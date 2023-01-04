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
 * @brief Explicitly computes the elements of the Schur complement represented
 * by `matrix`.
 * @param matrix A matrix created by `asSchurComplement`.
 */
template <class Policy>
Matrix<Policy> computeSchurComplement(const Matrix<Policy> &matrix);

extern template Matrix<SequentialComputePolicy>
computeSchurComplement(const Matrix<SequentialComputePolicy> &matrix);
extern template Matrix<ParallelComputePolicy>
computeSchurComplement(const Matrix<ParallelComputePolicy> &matrix);

} // namespace cpppetsc
} // namespace ae108

#include "ae108/cpppetsc/UniqueEntity.h"
#include <petscksp.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy>
Matrix<Policy> computeSchurComplement(const Matrix<Policy> &matrix) {
  auto result = Mat{};
  Policy::handleError(
      MatSchurComplementComputeExplicitOperator(matrix.data(), &result));
  return Matrix<Policy>(makeUniqueEntity<Policy>(result));
}

} // namespace cpppetsc
} // namespace ae108