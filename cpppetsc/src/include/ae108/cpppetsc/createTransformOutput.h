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
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Creates a vector that can store the result of a multiplication with
 * the provided matrix. It is initialized with zeroes.
 */
template <class Policy>
distributed<Vector<Policy>>
createTransformOutput(const Matrix<Policy> &transform) {
  auto out = Vec{};
  Policy::handleError(MatCreateVecs(transform.data(), nullptr, &out));
  auto result = distributed<Vector<Policy>>(makeUniqueEntity<Policy>(out));
  result.unwrap().setZero();
  return result;
}

extern template distributed<Vector<SequentialComputePolicy>>
createTransformOutput(const Matrix<SequentialComputePolicy> &);

extern template distributed<Vector<ParallelComputePolicy>>
createTransformOutput(const Matrix<ParallelComputePolicy> &);

} // namespace cpppetsc
} // namespace ae108