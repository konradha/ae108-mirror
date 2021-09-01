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

#include "ae108/cpppetsc/NonlinearSolver.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/cppptest/isLocal.h"
#include <gmock/gmock.h>

using ae108::cppptest::isLocal;
using ae108::cppptest::ValueNear;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct NonlinearSolver_Test : Test {
  using solver_type = NonlinearSolver<Policy>;
  using vector_type = typename solver_type::vector_type;
  using matrix_type = typename solver_type::matrix_type;

  /**
   * @brief Solves x^2 == 1 with the given solver.
   *
   * @return The solution.
   */
  distributed<vector_type>
  solveConvergingProblem(const solver_type &solver) const {
    return solver.solve(
        [](const distributed<vector_type> &input,
           distributed<vector_type> *const output) {
          const auto replacer = output->unwrap().replace();
          if (isLocal(output->unwrap(), 0)) {
            replacer(0) = std::norm(input(0)) - 1.;
          }
        },
        [](const distributed<vector_type> &input, matrix_type *const output) {
          const auto full = vector_type::fromDistributed(input);
          const auto replacer = output->assemblyView().replace();
          if (isLocal(*output, 0)) {
            replacer(0, 0) = 2. * full(0);
          }
        },
        tag<DistributedTag>(vector_type::fromList({7.})));
  }

  /**
   * @brief Solve f(x) = 1 with f'(x) = 1.
   */
  distributed<vector_type>
  solveNonConvergingProblem(const solver_type &solver) const {
    return solver.solve(
        [](const distributed<vector_type> &,
           distributed<vector_type> *const output) {
          const auto replacer = output->unwrap().replace();
          if (isLocal(output->unwrap(), 0)) {
            replacer(0) = 1.;
          }
        },
        [](const vector_type &, matrix_type *const output) {
          const auto replacer = output->assemblyView().replace();
          if (isLocal(*output, 0)) {
            replacer(0, 0) = 1.;
          }
        },
        tag<DistributedTag>(vector_type::fromList({7.})));
  }
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(NonlinearSolver_Test, Policies);

TYPED_TEST(NonlinearSolver_Test, solves_x2_eq_1) {
  using vector_type = typename TestFixture::solver_type::vector_type;

  typename TestFixture::solver_type solver(1);
  const auto solution =
      vector_type::fromDistributed(this->solveConvergingProblem(solver));

  EXPECT_THAT(solution(0), ValueNear(1., 1e-6));
}

TYPED_TEST(NonlinearSolver_Test, solves_x2_eq_1_when_provided_buffer) {
  using vector_type = typename TestFixture::solver_type::vector_type;
  using matrix_type = typename TestFixture::solver_type::matrix_type;

  typename TestFixture::solver_type solver(matrix_type(1, 1),
                                           tag<DistributedTag>(vector_type(1)));
  const auto solution =
      vector_type::fromDistributed(this->solveConvergingProblem(solver));

  EXPECT_THAT(solution(0), ValueNear(1., 1e-6));
}

TYPED_TEST(NonlinearSolver_Test,
           buffer_constructed_exits_with_exception_on_nonconvergence) {
  using vector_type = typename TestFixture::solver_type::vector_type;
  using matrix_type = typename TestFixture::solver_type::matrix_type;

  typename TestFixture::solver_type solver(matrix_type(1, 1),
                                           tag<DistributedTag>(vector_type(1)));
  EXPECT_THROW(this->solveNonConvergingProblem(solver),
               NonlinearSolverDivergedException);
}

TYPED_TEST(NonlinearSolver_Test,
           dimension_constructed_exits_with_exception_on_nonconvergence) {
  typename TestFixture::solver_type solver(1);
  EXPECT_THROW(this->solveNonConvergingProblem(solver),
               NonlinearSolverDivergedException);
}
} // namespace
} // namespace cpppetsc
} // namespace ae108
