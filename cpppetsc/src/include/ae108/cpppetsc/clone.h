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

#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/copy.h"

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Creates a clone of the vector.
 */
template <class Policy, class Tag>
TaggedEntity<Vector<Policy>, Tag>
clone(const TaggedEntity<Vector<Policy>, Tag> &vector) {
  auto duplicate = [&]() {
    auto vec = Vec{};
    Policy::handleError(VecDuplicate(vector.unwrap().data(), &vec));
    return TaggedEntity<Vector<Policy>, Tag>(makeUniqueEntity<Policy>(vec));
  }();
  copy(vector, &duplicate);
  return duplicate;
}

/**
 * @brief Creates a clone of the matrix.
 */
template <class Policy> Matrix<Policy> clone(const Matrix<Policy> &matrix) {
  auto mat = Mat{};
  Policy::handleError(MatDuplicate(matrix.data(), MAT_COPY_VALUES, &mat));
  return Matrix<Policy>(makeUniqueEntity<Policy>(mat));
}

extern template TaggedEntity<Vector<SequentialComputePolicy>, LocalTag>
clone(const TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> &);
extern template TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag>
clone(const TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> &);
extern template TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag>
clone(const TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> &);
extern template TaggedEntity<Vector<ParallelComputePolicy>, LocalTag>
clone(const TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> &);
extern template TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag>
clone(const TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> &);
extern template TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag>
clone(const TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> &);

extern template Matrix<SequentialComputePolicy>
clone(const Matrix<SequentialComputePolicy> &matrix);
extern template Matrix<ParallelComputePolicy>
clone(const Matrix<ParallelComputePolicy> &matrix);

} // namespace cpppetsc
} // namespace ae108