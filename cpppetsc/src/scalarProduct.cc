// © 2022 ETH Zurich, Mechanics and Materials Lab
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