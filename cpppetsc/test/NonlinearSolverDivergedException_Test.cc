// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#include "ae108/cpppetsc/NonlinearSolverDivergedException.h"
#include <gmock/gmock.h>

using testing::StrEq;
using testing::Test;

namespace ae108 {
namespace cpppetsc {
namespace {

struct NonlinearSolverDivergedException_Test : Test {
  NonlinearSolverDivergedException exception;
};

TEST_F(NonlinearSolverDivergedException_Test,
       what_states_that_nonlinear_solver_diverged) {
  EXPECT_THAT(exception.what(),
              StrEq("The nonlinear solver failed to find a solution."));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108
