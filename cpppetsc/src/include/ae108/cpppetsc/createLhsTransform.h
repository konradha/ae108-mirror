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
#include "ae108/cpppetsc/UniqueEntity.h"

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Creates a matrix P that can be multiplied from the left hand side
 * with the provided Matrix A, i.e. P * A is well-formed.
 */
template <class Policy>
Matrix<Policy>
createLhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::size_type rows) {
  using size_type = typename Matrix<Policy>::size_type;

  auto result = []() {
    auto mat = Mat{};
    Policy::handleError(MatCreate(Policy::communicator(), &mat));
    return Matrix<Policy>(makeUniqueEntity<Policy>(mat));
  }();

  const auto localCols = [&]() {
    auto rows = size_type{};
    auto cols = size_type{};
    Policy::handleError(MatGetLocalSize(matrix.data(), &rows, &cols));
    return rows;
  }();

  const auto columns = matrix.size().first;

  Policy::handleError(
      MatSetSizes(result.data(), PETSC_DECIDE, localCols, rows, columns));
  Policy::handleError(MatSetFromOptions(result.data()));
  Policy::handleError(MatSetUp(result.data()));
  result.finalize();
  return result;
}

extern template Matrix<SequentialComputePolicy>
createLhsTransform(const Matrix<SequentialComputePolicy> &matrix,
                   typename Matrix<SequentialComputePolicy>::size_type);
extern template Matrix<ParallelComputePolicy>
createLhsTransform(const Matrix<ParallelComputePolicy> &matrix,
                   typename Matrix<ParallelComputePolicy>::size_type);

} // namespace cpppetsc
} // namespace ae108