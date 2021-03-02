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
 * @brief Creates a vector that can be multiplied with the provided matrix.
 * It is initialized with zeroes.
 */
template <class Policy>
distributed<Vector<Policy>>
createTransformInput(const Matrix<Policy> &transform) {
  auto in = Vec{};
  Policy::handleError(MatCreateVecs(transform.data(), &in, nullptr));
  auto result = distributed<Vector<Policy>>(makeUniqueEntity<Policy>(in));
  result.unwrap().setZero();
  return result;
}

extern template distributed<Vector<SequentialComputePolicy>>
createTransformInput(const Matrix<SequentialComputePolicy> &);

extern template distributed<Vector<ParallelComputePolicy>>
createTransformInput(const Matrix<ParallelComputePolicy> &);

} // namespace cpppetsc
} // namespace ae108