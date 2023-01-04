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