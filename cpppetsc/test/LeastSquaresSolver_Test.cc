// © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/LeastSquaresSolver.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/cppptest/isLocal.h"
#include <gmock/gmock.h>

using ae108::cppptest::isLocal;
using ae108::cppptest::ScalarNear;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct LeastSquaresSolver_Test : Test {
  using solver_type = LeastSquaresSolver<Policy>;
  using vector_type = typename solver_type::vector_type;
  using matrix_type = typename solver_type::matrix_type;

  /**
   * @brief Solves {x == 1, 2x == 4} by minimizing the squared distance.
   *
   * @return The solution.
   */
  distributed<vector_type>
  solveConflictingEquations(const solver_type &solver) const {
    return solver.solve(
        [](const distributed<vector_type> &input,
           distributed<vector_type> *const output) {
          const auto full = vector_type::fromDistributed(input);
          const auto replacer = output->unwrap().replace();
          if (isLocal(output->unwrap(), 0)) {
            replacer(0) = full(0) - 1.;
          }
          if (isLocal(output->unwrap(), 1)) {
            replacer(1) = 2. * full(0) - 4.;
          }
        },
        [](const distributed<vector_type> &, matrix_type *const output) {
          const auto replacer = output->assemblyView().replace();
          if (isLocal(*output, 0)) {
            replacer(0, 0) = 1.;
          }
          if (isLocal(*output, 1)) {
            replacer(1, 0) = 2.;
          }
        },
        tag<DistributedTag>(vector_type::fromList({7.})));
  }
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(LeastSquaresSolver_Test, Policies);

TYPED_TEST(LeastSquaresSolver_Test,
           finds_optimal_solution_for_conflicting_equations) {
  using matrix_type = typename TestFixture::solver_type::matrix_type;
  using vector_type = typename TestFixture::solver_type::vector_type;

  typename TestFixture::solver_type solver(matrix_type(2, 1),
                                           distributed<vector_type>(2));
  const auto solution =
      vector_type::fromDistributed(this->solveConflictingEquations(solver));

  EXPECT_THAT(solution(0), ScalarNear(18. / 10., 1e-6));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108
