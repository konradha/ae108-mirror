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