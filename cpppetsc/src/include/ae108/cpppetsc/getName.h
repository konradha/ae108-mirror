// © 2020 ETH Zurich, Mechanics and Materials Lab
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

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Returns the name of the provided vector.
 */
template <class Policy>
const char *getName(const distributed<Vector<Policy>> &vector) {
  const char *result = nullptr;
  Policy::handleError(
      PetscObjectGetName((PetscObject)vector.unwrap().data(), &result));
  return result;
}

extern template const char *getName<SequentialComputePolicy>(
    const distributed<Vector<SequentialComputePolicy>> &);
extern template const char *getName<ParallelComputePolicy>(
    const distributed<Vector<ParallelComputePolicy>> &);

} // namespace cpppetsc
} // namespace ae108