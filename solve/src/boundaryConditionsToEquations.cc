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

#include "ae108/solve/boundaryConditionsToEquations.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace solve {

template AffineTransform<cpppetsc::SequentialComputePolicy>
boundaryConditionsToEquations(
    const std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<
        cpppetsc::Mesh<cpppetsc::SequentialComputePolicy>>> &,
    const cpppetsc::Mesh<cpppetsc::SequentialComputePolicy> &);

template AffineTransform<cpppetsc::ParallelComputePolicy>
boundaryConditionsToEquations(
    const std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<
        cpppetsc::Mesh<cpppetsc::ParallelComputePolicy>>> &,
    const cpppetsc::Mesh<cpppetsc::ParallelComputePolicy> &);

} // namespace solve
} // namespace ae108