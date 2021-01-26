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

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Copies the contents of vector `from` to vector `to`. The layouts of
 * the vectors must be identical.
 */
template <class Policy, class Tag>
void copy(const TaggedEntity<Vector<Policy>, Tag> &from,
          TaggedEntity<Vector<Policy>, Tag> *to) {
  Policy::handleError(VecCopy(from.unwrap().data(), to->unwrap().data()));
}

extern template void
copy(const TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> &,
     TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> *);
extern template void
copy(const TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> &,
     TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> *);
extern template void
copy(const TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> &,
     TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> *);

extern template void
copy(const TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> &,
     TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> *);
extern template void
copy(const TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> &,
     TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> *);
extern template void
copy(const TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> &,
     TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> *);

} // namespace cpppetsc
} // namespace ae108