// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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
