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

#include "ae108/cpppetsc/copy.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace cpppetsc {

template void
copy(const TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> &,
     TaggedEntity<Vector<SequentialComputePolicy>, LocalTag> *);
template void
copy(const TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> &,
     TaggedEntity<Vector<SequentialComputePolicy>, DistributedTag> *);
template void
copy(const TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> &,
     TaggedEntity<Vector<SequentialComputePolicy>, GlobalTag> *);

template void
copy(const TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> &,
     TaggedEntity<Vector<ParallelComputePolicy>, LocalTag> *);
template void
copy(const TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> &,
     TaggedEntity<Vector<ParallelComputePolicy>, DistributedTag> *);
template void
copy(const TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> &,
     TaggedEntity<Vector<ParallelComputePolicy>, GlobalTag> *);

} // namespace cpppetsc
} // namespace ae108