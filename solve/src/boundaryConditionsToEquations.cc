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