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

#include "ae108/cpppetsc/clone.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace cpppetsc {

template TaggedEntity<Vector<SequentialComputePolicy>, LocalTag>
clone(const TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> &);
template TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag>
clone(const TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> &);
template TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag>
clone(const TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> &);
template TaggedEntity<Vector<ParallelComputePolicy>, LocalTag>
clone(const TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> &);
template TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag>
clone(const TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> &);
template TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag>
clone(const TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> &);

template Matrix<SequentialComputePolicy>
clone(const Matrix<SequentialComputePolicy> &matrix);
template Matrix<ParallelComputePolicy>
clone(const Matrix<ParallelComputePolicy> &matrix);

} // namespace cpppetsc
} // namespace ae108