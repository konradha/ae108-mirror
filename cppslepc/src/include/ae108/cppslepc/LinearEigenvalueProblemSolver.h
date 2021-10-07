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

#pragma once

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include <slepceps.h>

namespace ae108 {
namespace cppslepc {

template <class Policy> class LinearEigenvalueProblemSolver {
public:
  using size_type = PetscInt;
  using value_type = PetscScalar;
  using real_type = PetscReal;

  explicit LinearEigenvalueProblemSolver();

  /**
   * @brief Returns the internal solver.
   */
  EPS data() const;

private:
  cpppetsc::UniqueEntity<EPS> eps_;
};

extern template class LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>;
extern template class LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>;

} // namespace cppslepc
} // namespace ae108

namespace ae108 {
namespace cppslepc {

template <class Policy>
LinearEigenvalueProblemSolver<Policy>::LinearEigenvalueProblemSolver()
    : eps_([]() {
        auto solver = EPS{};
        Policy::handleError(EPSCreate(Policy::communicator(), &solver));
        return cpppetsc::UniqueEntity<EPS>(
            solver, [](EPS eps) { Policy::handleError(EPSDestroy(&eps)); });
      }()) {
  Policy::handleError(EPSSetFromOptions(eps_.get()));
}

template <class Policy>
EPS LinearEigenvalueProblemSolver<Policy>::data() const {
  return eps_.get();
}

} // namespace cppslepc
} // namespace ae108