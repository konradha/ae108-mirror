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

#pragma once

#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Computes the scalar product of the two vectors `x` and `y`.
 */
template <class Policy, class Tag>
typename Vector<Policy>::value_type
scalarProduct(const TaggedEntity<Vector<Policy>, Tag> &x,
              const TaggedEntity<Vector<Policy>, Tag> &y) {
  auto result = typename Vector<Policy>::value_type();
  Policy::handleError(VecDot(x.unwrap().data(), y.unwrap().data(), &result));
  return result;
}

} // namespace cpppetsc
} // namespace ae108