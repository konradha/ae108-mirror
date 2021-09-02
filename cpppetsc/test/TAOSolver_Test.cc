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

#ifndef AE108_PETSC_COMPLEX

#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TAOSolver.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/cppptest/isLocal.h"
#include <cmath>
#include <gmock/gmock.h>

using ae108::cppptest::AlmostEqIfLocal;
using ae108::cppptest::isLocal;
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
  using value_type = typename TestFixture::value_type;
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

  EXPECT_THAT(solution.unwrap(), AlmostEqIfLocal(0, 1.));
}

TYPED_TEST(TAOSolver_Test, minimizes_two_dimensional_problem) {
  using value_type = typename TestFixture::value_type;
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

  EXPECT_THAT(solution.unwrap(), AlmostEqIfLocal(0, 1.));
  EXPECT_THAT(solution.unwrap(), AlmostEqIfLocal(1, 2.));
}

TYPED_TEST(TAOSolver_Test, raises_exception_on_nonconvergence) {
  using value_type = typename TestFixture::value_type;
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
  EXPECT_THAT(solution.unwrap(), AlmostEqIfLocal(0, 2.));
}

TYPED_TEST(TAOSolver_BoundsTest, minimizes_with_upper_bounds) {
  using vector_type = typename TestFixture::vector_type;
  this->solver.setBounds(tag<DistributedTag>(vector_type::fromList(
                             {TestFixture::solver_type::no_lower_bound})),
                         tag<DistributedTag>(vector_type::fromList({-2.})));

  const auto solution = this->minimizeQuadratic(-7.);
  EXPECT_THAT(solution.unwrap(), AlmostEqIfLocal(0, -2.));
}

TYPED_TEST(TAOSolver_BoundsTest, minimizes_with_both_bounds) {
  using vector_type = typename TestFixture::vector_type;
  this->solver.setBounds(tag<DistributedTag>(vector_type::fromList({2.})),
                         tag<DistributedTag>(vector_type::fromList({2.})));

  const auto solution = this->minimizeQuadratic(-7.);
  EXPECT_THAT(solution.unwrap(), AlmostEqIfLocal(0, 2.));
}

TYPED_TEST(TAOSolver_BoundsTest, lmvm_does_not_apply_bounds) {
  this->solver.setType(TestFixture::solver_type::Type::lmvm);
  using vector_type = typename TestFixture::vector_type;
  this->solver.setBounds(tag<DistributedTag>(vector_type::fromList({2.})),
                         tag<DistributedTag>(vector_type::fromList({2.})));

  const auto solution = this->minimizeQuadratic(-7.);
  EXPECT_THAT(solution.unwrap(), AlmostEqIfLocal(0, 1.));
}
} // namespace
} // namespace cpppetsc
} // namespace ae108

#endif