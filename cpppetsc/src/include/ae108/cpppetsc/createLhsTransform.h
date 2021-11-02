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
 *
 * @param rows Global number of rows of the created matrix.
 */
template <class Policy>
Matrix<Policy>
createLhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::size_type rows);

/**
 * @brief Creates a matrix P that can be multiplied from the left hand side
 * with the provided Matrix A, i.e. P * A is well-formed.
 */
template <class Policy>
Matrix<Policy>
createLhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::LocalRows localRows,
                   const typename Matrix<Policy>::GlobalRows globalRows);

extern template Matrix<SequentialComputePolicy>
createLhsTransform(const Matrix<SequentialComputePolicy> &,
                   typename Matrix<SequentialComputePolicy>::size_type);
extern template Matrix<ParallelComputePolicy>
createLhsTransform(const Matrix<ParallelComputePolicy> &,
                   typename Matrix<ParallelComputePolicy>::size_type);

extern template Matrix<SequentialComputePolicy>
createLhsTransform(const Matrix<SequentialComputePolicy> &,
                   const typename Matrix<SequentialComputePolicy>::LocalRows,
                   const typename Matrix<SequentialComputePolicy>::GlobalRows);
extern template Matrix<ParallelComputePolicy>
createLhsTransform(const Matrix<ParallelComputePolicy> &,
                   const typename Matrix<ParallelComputePolicy>::LocalRows,
                   const typename Matrix<ParallelComputePolicy>::GlobalRows);

} // namespace cpppetsc
} // namespace ae108

namespace ae108 {
namespace cpppetsc {

template <class Policy>
Matrix<Policy>
createLhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::size_type rows) {
  return createLhsTransform(matrix,
                            typename Matrix<Policy>::LocalRows{PETSC_DECIDE},
                            typename Matrix<Policy>::GlobalRows{rows});
}

template <class Policy>
Matrix<Policy>
createLhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::LocalRows localRows,
                   const typename Matrix<Policy>::GlobalRows globalRows) {
  return Matrix<Policy>(
      localRows, typename Matrix<Policy>::LocalCols{matrix.localSize().first},
      globalRows, typename Matrix<Policy>::GlobalCols{matrix.size().first});
}

} // namespace cpppetsc
} // namespace ae108