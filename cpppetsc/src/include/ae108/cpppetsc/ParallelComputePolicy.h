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

#pragma once

#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include <mpi.h>
#include <petscsys.h>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief The policies configure the wrapper classes (error handling,
 * communicator).
 */
struct ParallelComputePolicy {
  static inline MPI_Comm communicator();

  static inline void handleError(const PetscErrorCode code);

  static inline bool isPrimaryRank();
};

/********************************************************************
 *  implementations
 *******************************************************************/

MPI_Comm ParallelComputePolicy::communicator() { return PETSC_COMM_WORLD; }

void ParallelComputePolicy::handleError(const PetscErrorCode code) {
  CHKERRABORT(communicator(), code);
}

bool ParallelComputePolicy::isPrimaryRank() {
  int rank = 0;
  handleError(MPI_Comm_rank(communicator(), &rank));
  return rank == 0;
}
} // namespace cpppetsc
} // namespace ae108