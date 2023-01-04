// © 2020 ETH Zurich, Mechanics and Materials Lab
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
 * @brief Returns the name of the provided vector.
 */
template <class Policy>
const char *getName(const distributed<Vector<Policy>> &vector) {
  const char *result = nullptr;
  Policy::handleError(
      PetscObjectGetName((PetscObject)vector.unwrap().data(), &result));
  return result;
}

extern template const char *getName<SequentialComputePolicy>(
    const distributed<Vector<SequentialComputePolicy>> &);
extern template const char *getName<ParallelComputePolicy>(
    const distributed<Vector<ParallelComputePolicy>> &);

} // namespace cpppetsc
} // namespace ae108