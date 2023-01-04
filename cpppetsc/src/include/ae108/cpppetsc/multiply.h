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
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Multiplies the given matrix by the given vector and returns the
 * result.
 */
template <class Policy>
distributed<Vector<Policy>> multiply(const Matrix<Policy> &matrix,
                                     const distributed<Vector<Policy>> &vector);

extern template distributed<Vector<SequentialComputePolicy>>
multiply(const Matrix<SequentialComputePolicy> &matrix,
         const distributed<Vector<SequentialComputePolicy>> &vector);

extern template distributed<Vector<ParallelComputePolicy>>
multiply(const Matrix<ParallelComputePolicy> &matrix,
         const distributed<Vector<ParallelComputePolicy>> &vector);

} // namespace cpppetsc
} // namespace ae108

#include "ae108/cpppetsc/createTransformOutput.h"
#include <petscmat.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy>
distributed<Vector<Policy>>
multiply(const Matrix<Policy> &matrix,
         const distributed<Vector<Policy>> &vector) {
  auto output = createTransformOutput(matrix);
  Policy::handleError(
      MatMult(matrix.data(), vector.unwrap().data(), output.unwrap().data()));
  return output;
}

} // namespace cpppetsc
} // namespace ae108