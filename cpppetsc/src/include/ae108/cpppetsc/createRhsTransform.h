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
 * @brief Creates a matrix P that can be multiplied from the right hand side
 * with the provided Matrix A, i.e. A * P is well-formed.
 *
 * @param columns Global number of columns of the created matrix.
 */
template <class Policy>
Matrix<Policy>
createRhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::size_type columns);

/**
 * @brief Creates a matrix P that can be multiplied from the right hand side
 * with the provided Matrix A, i.e. A * P is well-formed.
 */
template <class Policy>
Matrix<Policy>
createRhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::LocalCols localColumns,
                   const typename Matrix<Policy>::GlobalCols globalColumns);

extern template Matrix<SequentialComputePolicy>
createRhsTransform(const Matrix<SequentialComputePolicy> &,
                   typename Matrix<SequentialComputePolicy>::size_type);
extern template Matrix<ParallelComputePolicy>
createRhsTransform(const Matrix<ParallelComputePolicy> &,
                   typename Matrix<ParallelComputePolicy>::size_type);

extern template Matrix<SequentialComputePolicy>
createRhsTransform(const Matrix<SequentialComputePolicy> &,
                   const typename Matrix<SequentialComputePolicy>::LocalCols,
                   const typename Matrix<SequentialComputePolicy>::GlobalCols);
extern template Matrix<ParallelComputePolicy>
createRhsTransform(const Matrix<ParallelComputePolicy> &,
                   const typename Matrix<ParallelComputePolicy>::LocalCols,
                   const typename Matrix<ParallelComputePolicy>::GlobalCols);

} // namespace cpppetsc
} // namespace ae108

namespace ae108 {
namespace cpppetsc {

template <class Policy>
Matrix<Policy>
createRhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::size_type columns) {
  return createRhsTransform(matrix,
                            typename Matrix<Policy>::LocalCols{PETSC_DECIDE},
                            typename Matrix<Policy>::GlobalCols{columns});
}

template <class Policy>
Matrix<Policy>
createRhsTransform(const Matrix<Policy> &matrix,
                   const typename Matrix<Policy>::LocalCols localColumns,
                   const typename Matrix<Policy>::GlobalCols globalColumns) {
  return Matrix<Policy>(
      typename Matrix<Policy>::LocalRows{matrix.localSize().second},
      localColumns, typename Matrix<Policy>::GlobalRows{matrix.size().second},
      globalColumns);
}

} // namespace cpppetsc
} // namespace ae108