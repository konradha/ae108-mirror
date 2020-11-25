// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#include "ae108/cpppetsc/createVectorFromSource.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace cpppetsc {

template distributed<typename Mesh<SequentialComputePolicy>::vector_type>
createVectorFromSource(
    const Mesh<SequentialComputePolicy> &,
    typename Mesh<SequentialComputePolicy>::size_type,
    std::function<void(typename Mesh<SequentialComputePolicy>::size_type,
                       typename Mesh<SequentialComputePolicy>::value_type *)>);

template distributed<typename Mesh<ParallelComputePolicy>::vector_type>
createVectorFromSource(
    const Mesh<ParallelComputePolicy> &,
    typename Mesh<ParallelComputePolicy>::size_type,
    std::function<void(typename Mesh<ParallelComputePolicy>::size_type,
                       typename Mesh<ParallelComputePolicy>::value_type *)>);

} // namespace cpppetsc
} // namespace ae108
