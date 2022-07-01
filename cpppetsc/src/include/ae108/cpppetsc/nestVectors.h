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