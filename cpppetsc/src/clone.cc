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