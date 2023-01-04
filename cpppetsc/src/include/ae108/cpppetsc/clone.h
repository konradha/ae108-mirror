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

#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/copy.h"
#include <petscmat.h>
#include <petscvec.h>

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