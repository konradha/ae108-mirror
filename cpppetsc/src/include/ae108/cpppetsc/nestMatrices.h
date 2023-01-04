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