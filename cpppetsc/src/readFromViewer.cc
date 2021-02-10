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

#include "ae108/cpppetsc/readFromViewer.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace cpppetsc {

template void
readFromViewer(const Viewer<SequentialComputePolicy> &viewer,
               distributed<Vector<SequentialComputePolicy>> *const vector);
template void
readFromViewer(const Viewer<ParallelComputePolicy> &viewer,
               distributed<Vector<ParallelComputePolicy>> *const vector);

} // namespace cpppetsc
} // namespace ae108