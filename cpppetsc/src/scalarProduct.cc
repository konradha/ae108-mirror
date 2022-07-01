// © 2022 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/scalarProduct.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace cpppetsc {

template typename Vector<SequentialComputePolicy>::value_type
scalarProduct(const TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> &,
              const TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> &);
template typename Vector<SequentialComputePolicy>::value_type scalarProduct(
    const TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> &,
    const TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> &);
template typename Vector<SequentialComputePolicy>::value_type
scalarProduct(const TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> &,
              const TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> &);

template typename Vector<ParallelComputePolicy>::value_type
scalarProduct(const TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> &,
              const TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> &);
template typename Vector<ParallelComputePolicy>::value_type scalarProduct(
    const TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> &,
    const TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> &);
template typename Vector<ParallelComputePolicy>::value_type
scalarProduct(const TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> &,
              const TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> &);
} // namespace cpppetsc
} // namespace ae108