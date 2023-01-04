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

#include "ae108/cpppetsc/LinearSolver.h"
#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct LinearSolver_Test : Test {
  using linearsolver_type = LinearSolver<Policy>;
  using vector_type = typename linearsolver_type::vector_type;
  using matrix_type = typename linearsolver_type::matrix_type;

  matrix_type matrix = matrix_type::fromList({{2., 0.}, {0., 3.}});
  distributed<vector_type> vector =
      tag<DistributedTag>(vector_type::fromList({8., 6.}));

  linearsolver_type solver = linearsolver_type(&matrix);
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(LinearSolver_Test, Policies);

TYPED_TEST(LinearSolver_Test, solving_works) {
  const auto result = TestFixture::vector_type::fromDistributed(
      this->solver.solve(this->vector));

  EXPECT_THAT(result(0), ScalarEq(4.));
  EXPECT_THAT(result(1), ScalarEq(2.));
}

TYPED_TEST(LinearSolver_Test, solving_works_with_provided_result_vector) {
  const auto result =
      TestFixture::vector_type::fromDistributed(this->solver.solve(
          this->vector,
          tag<DistributedTag>(
              TestFixture::vector_type::fromLayoutOf(this->vector))));

  EXPECT_THAT(result(0), ScalarEq(4.));
  EXPECT_THAT(result(1), ScalarEq(2.));
}

TYPED_TEST(LinearSolver_Test, move_construction_works) {
  const auto copy(std::move(this->solver));

  const auto result =
      TestFixture::vector_type::fromDistributed(copy.solve(this->vector));

  EXPECT_THAT(result(0), ScalarEq(4.));
  EXPECT_THAT(result(1), ScalarEq(2.));
}

TYPED_TEST(LinearSolver_Test, move_assignment_works) {
  const auto matrix_2 =
      TestFixture::matrix_type::fromList({{1., 0.}, {0., 1.}});
  typename TestFixture::linearsolver_type copy(&matrix_2);

  copy = std::move(this->solver);
  const auto result =
      TestFixture::vector_type::fromDistributed(copy.solve(this->vector));

  EXPECT_THAT(result(0), ScalarEq(4.));
  EXPECT_THAT(result(1), ScalarEq(2.));
}

TYPED_TEST(LinearSolver_Test, solving_singular_system_reports_error) {
  const auto matrix_2 =
      TestFixture::matrix_type::fromList({{0., 0.}, {0., 0.}});
  const typename TestFixture::linearsolver_type solver(&matrix_2);

  EXPECT_THROW(solver.solve(this->vector), LinearSolverDivergedException);
}

TYPED_TEST(LinearSolver_Test, solver_keeps_reference) {
  auto matrix_2 = TestFixture::matrix_type::fromList({{0., 0.}, {0., 0.}});
  const typename TestFixture::linearsolver_type solver(&matrix_2);

  {
    auto replacer = matrix_2.assemblyView().replace();
    const auto range = matrix_2.localRowRange();
    if (range.first <= 0 && 0 < range.second) {
      replacer.element(0, 0, 1.);
    }
    if (range.first <= 1 && 1 < range.second) {
      replacer.element(1, 1, 1.);
    }
  }

  const auto result =
      TestFixture::vector_type::fromDistributed(solver.solve(this->vector));

  EXPECT_THAT(result(0), ScalarEq(8.));
  EXPECT_THAT(result(1), ScalarEq(6.));
}
} // namespace
} // namespace cpppetsc
} // namespace ae108
