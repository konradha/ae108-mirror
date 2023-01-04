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