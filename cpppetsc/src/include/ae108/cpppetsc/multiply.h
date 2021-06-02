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