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