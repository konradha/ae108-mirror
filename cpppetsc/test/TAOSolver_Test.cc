// © 2020, 2022 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TAOSolver.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/cppptest/isLocal.h"
#include <cmath>
#include <gmock/gmock.h>

using ae108::cppptest::isLocal;
using ae108::cppptest::ScalarEqIfLocal;
using ae108::cppptest::ScalarNearIfLocal;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct TAOSolver_Test : Test {
  using solver_type = TAOSolver<Policy>;
  using value_type = typename solver_type::value_type;
  using real_type = typename solver_type::real_type;
  using vector_type = typename solver_type::vector_type;
  using matrix_type = typename solver_type::matrix_type;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(TAOSolver_Test, Policies);

TYPED_TEST(TAOSolver_Test, minimizes_x_minus_1_squared) {
  using real_type = typename TestFixture::real_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;

  typename TestFixture::solver_type solver(matrix_type(1, 1));
  const auto solution = solver.solve(
      [](const distributed<vector_type> &input, real_type *const output) {
        const auto full = vector_type::fromDistributed(input);
        *output = std::norm(full(0) - 1.);
      },
      [](const distributed<vector_type> &input,
         distributed<vector_type> *const output) {
        const auto full = vector_type::fromDistributed(input);
        const auto replacer = output->unwrap().replace();
        if (isLocal(output->unwrap(), 0)) {
          replacer(0) = 2. * (full(0) - 1.);
        }
      },
      [](const distributed<vector_type> &, matrix_type *const output) {
        const auto replacer = output->assemblyView().replace();
        if (isLocal(*output, 0)) {
          replacer(0, 0) = 2.;
        }
      },
      tag<DistributedTag>(vector_type::fromList({7.})));

  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 1.));
}

TYPED_TEST(TAOSolver_Test, minimizes_two_dimensional_problem) {
  using real_type = typename TestFixture::real_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;

  typename TestFixture::solver_type solver(matrix_type(2, 2));
  const auto solution = solver.solve(
      [](const distributed<vector_type> &input, real_type *const output) {
        const auto full = vector_type::fromDistributed(input);
        *output =
            std::norm(full(0) - 1) + std::norm(2. * full(0) - full(1)) + 7;
      },
      [](const distributed<vector_type> &input,
         distributed<vector_type> *const output) {
        const auto full = vector_type::fromDistributed(input);
        const auto replacer = output->unwrap().replace();
        if (isLocal(output->unwrap(), 0)) {
          replacer(0) = 2. * (full(0) - 1.) + 4. * (2. * full(0) - full(1));
        }
        if (isLocal(output->unwrap(), 1)) {
          replacer(1) = -2. * (2. * full(0) - full(1));
        }
      },
      [](const distributed<vector_type> &, matrix_type *const output) {
        const auto replacer = output->assemblyView().replace();
        if (isLocal(*output, 0)) {
          replacer(0, 0) = 10.;
          replacer(0, 1) = -4.;
        }
        if (isLocal(*output, 1)) {
          replacer(1, 0) = -4.;
          replacer(1, 1) = 2.;
        }
      },
      tag<DistributedTag>(vector_type::fromList({3., 4.})));

  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 1.));
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(1, 2.));
}

TYPED_TEST(TAOSolver_Test, raises_exception_on_nonconvergence) {
  using real_type = typename TestFixture::real_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;

  typename TestFixture::solver_type solver(matrix_type(1, 1));
  const auto callSolve = [&solver]() {
    solver.solve(
        [](const distributed<vector_type> &, real_type *const output) {
          *output = 1.;
        },
        [](const distributed<vector_type> &,
           distributed<vector_type> *const output) {
          const auto replacer = output->unwrap().replace();
          replacer(0) = 1.;
        },
        [](const distributed<vector_type> &, matrix_type *const output) {
          const auto replacer = output->assemblyView().replace();
          replacer(0, 0) = 0.;
        },
        tag<DistributedTag>(vector_type::fromList({1.})));
  };

  EXPECT_THROW(callSolve(), TAOSolverDivergedException);
}

TYPED_TEST(TAOSolver_Test, minimizes_with_equality_constraints) {
  using real_type = typename TestFixture::real_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;

  typename TestFixture::solver_type solver(
      matrix_type(2, 2), TestFixture::solver_type::Type::almm);

  const auto value = 5.;

  solver.setConstraints({
      [&](const distributed<vector_type> &input,
          distributed<vector_type> *const output) {
        const auto full = vector_type::fromDistributed(input);
        output->unwrap().setZero();
        output->unwrap().replace().element(0, full(0) - value);
      },
      [](const distributed<vector_type> &, matrix_type *const output) {
        output->setZero();
        const auto replacer =
            output->assemblyView().replace().element(0, 0, 1.);
      },
      {tag<DistributedTag>(vector_type(1)), matrix_type(1, 2)},
  });

  const auto solution = solver.solve(
      [](const distributed<vector_type> &input, real_type *const output) {
        const auto full = vector_type::fromDistributed(input);
        *output = std::pow(full(0) - 1., 2.) + std::pow(full(1) - 2., 2.);
      },
      [](const distributed<vector_type> &input,
         distributed<vector_type> *const output) {
        const auto full = vector_type::fromDistributed(input);
        const auto replacer = output->unwrap().replace();
        replacer(0) = 2. * (full(0) - 1.);
        replacer(1) = 2. * (full(1) - 2.);
      },
      [](const distributed<vector_type> &, matrix_type *const output) {
        const auto replacer = output->assemblyView().replace();
        replacer(0, 0) = 2.;
        replacer(1, 1) = 2.;
      },
      tag<DistributedTag>(vector_type::fromList({3., 4.})));

  EXPECT_THAT(solution.unwrap(), ScalarNearIfLocal(0, value, 1e-7));
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(1, 2.));
}

template <class Policy> struct TAOSolver_BoundsTest : Test {
  using solver_type = TAOSolver<Policy>;
  using value_type = typename solver_type::value_type;
  using real_type = typename solver_type::real_type;
  using vector_type = typename solver_type::vector_type;
  using matrix_type = typename solver_type::matrix_type;

  solver_type solver{matrix_type(1, 1), solver_type::Type::blmvm};

  distributed<vector_type> minimizeQuadratic(const value_type guess) {
    return solver.solve(
        [](const distributed<vector_type> &input, real_type *const output) {
          const auto full = vector_type::fromDistributed(input);
          *output = std::norm(full(0) - 1.);
        },
        [](const distributed<vector_type> &input,
           distributed<vector_type> *const output) {
          const auto full = vector_type::fromDistributed(input);
          const auto replacer = output->unwrap().replace();
          if (isLocal(output->unwrap(), 0)) {
            replacer(0) = 2. * (full(0) - 1.);
          }
        },
        [](const distributed<vector_type> &, matrix_type *const output) {
          const auto replacer = output->assemblyView().replace();
          if (isLocal(*output, 0)) {
            replacer(0, 0) = 2.;
          }
        },
        tag<DistributedTag>(vector_type::fromList({guess})));
  }
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(TAOSolver_BoundsTest, Policies);

TYPED_TEST(TAOSolver_BoundsTest, minimizes_with_lower_bounds) {
  using vector_type = typename TestFixture::vector_type;
  this->solver.setBounds(tag<DistributedTag>(vector_type::fromList({2.})),
                         tag<DistributedTag>(vector_type::fromList(
                             {TestFixture::solver_type::no_upper_bound})));

  const auto solution = this->minimizeQuadratic(7.);
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 2.));
}

TYPED_TEST(TAOSolver_BoundsTest, minimizes_with_upper_bounds) {
  using vector_type = typename TestFixture::vector_type;
  this->solver.setBounds(tag<DistributedTag>(vector_type::fromList(
                             {TestFixture::solver_type::no_lower_bound})),
                         tag<DistributedTag>(vector_type::fromList({-2.})));

  const auto solution = this->minimizeQuadratic(-7.);
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, -2.));
}

TYPED_TEST(TAOSolver_BoundsTest, minimizes_with_both_bounds) {
  using vector_type = typename TestFixture::vector_type;
  this->solver.setBounds(tag<DistributedTag>(vector_type::fromList({2.})),
                         tag<DistributedTag>(vector_type::fromList({2.})));

  const auto solution = this->minimizeQuadratic(-7.);
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 2.));
}

TYPED_TEST(TAOSolver_BoundsTest, lmvm_does_not_apply_bounds) {
  this->solver.setType(TestFixture::solver_type::Type::lmvm);
  using vector_type = typename TestFixture::vector_type;
  this->solver.setBounds(tag<DistributedTag>(vector_type::fromList({2.})),
                         tag<DistributedTag>(vector_type::fromList({2.})));

  const auto solution = this->minimizeQuadratic(-7.);
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 1.));
}
} // namespace
} // namespace cpppetsc
} // namespace ae108

#endif