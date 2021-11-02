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

#include "ae108/cpppetsc/InvalidParametersException.h"
#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include <algorithm>
#include <vector>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Creates a matrix that consists of all the provided matrices combined.
 * The matrices are not actually copied; instead the returned matrix contains a
 * reference to the provided matrices.
 *
 * @throw InvalidParametersException if the number of columns is inconsistent.
 */
template <class Policy>
Matrix<Policy>
nestMatrices(const std::vector<std::vector<Matrix<Policy> *>> &matrices) {
  using matrix_type = Matrix<Policy>;

  const auto rows = matrices.size();
  const auto cols = matrices.empty() ? 0 : matrices.front().size();

  std::vector<Mat> data;
  data.reserve(rows * cols);
  for (auto &&row : matrices) {
    if (row.size() != cols) {
      throw InvalidParametersException{};
    }
    for (auto &&col : row) {
      data.push_back(col->data());
    }
  }

  auto nested = Mat{};
  Policy::handleError(MatCreateNest(Policy::communicator(), rows, nullptr, cols,
                                    nullptr, data.data(), &nested));
  return matrix_type(cpppetsc::makeUniqueEntity<Policy>(nested));
}

extern template Matrix<SequentialComputePolicy> nestMatrices(
    const std::vector<std::vector<Matrix<SequentialComputePolicy> *>> &);

extern template Matrix<ParallelComputePolicy>
nestMatrices(const std::vector<std::vector<Matrix<ParallelComputePolicy> *>> &);

} // namespace cpppetsc
} // namespace ae108