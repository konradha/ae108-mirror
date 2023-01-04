// © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
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

#ifndef AE108_PETSC_COMPLEX

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

#endif