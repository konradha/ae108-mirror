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
#include <algorithm>
#include <vector>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Creates a vector that consists of all the provided vectors stacked on
 * top of each other. The vectors are not actually copied; instead the returned
 * vector contains a reference to the provided vectors.
 */
template <class Policy>
distributed<Vector<Policy>> nestVectors(
    const std::vector<cpppetsc::distributed<Vector<Policy>> *> &vectors) {
  using vector_type = Vector<Policy>;

  std::vector<Vec> data(vectors.size());
  std::transform(
      vectors.begin(), vectors.end(), data.begin(),
      [](cpppetsc::distributed<vector_type> *v) { return v->unwrap().data(); });

  auto nested = Vec{};
  Policy::handleError(
      VecCreateNest(Policy::communicator(),
                    static_cast<typename vector_type::size_type>(data.size()),
                    nullptr, data.data(), &nested));
  return cpppetsc::distributed<vector_type>(
      cpppetsc::makeUniqueEntity<Policy>(nested));
}

extern template distributed<Vector<SequentialComputePolicy>> nestVectors(
    const std::vector<cpppetsc::distributed<Vector<SequentialComputePolicy>> *>
        &);

extern template distributed<Vector<ParallelComputePolicy>> nestVectors(
    const std::vector<cpppetsc::distributed<Vector<ParallelComputePolicy>> *>
        &);

} // namespace cpppetsc
} // namespace ae108