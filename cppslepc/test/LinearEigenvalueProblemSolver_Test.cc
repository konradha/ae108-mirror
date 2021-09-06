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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cppslepc/LinearEigenvalueProblemSolver.h"
#include <gmock/gmock.h>

using testing::ElementsAre;
using testing::Eq;
using testing::Test;
using testing::Types;

namespace {
MATCHER_P2(ComplexNear, reference, tolerance,
           (negation ? std::string("is") : std::string("isn't")) +
               " equal to " + ::testing::PrintToString(reference)) {
  return ::testing::ExplainMatchResult(
             ::testing::DoubleNear(std::complex<double>(reference).real(),
                                   tolerance),
             std::complex<double>(arg).real(), result_listener) &&
         ::testing::ExplainMatchResult(
             ::testing::DoubleNear(std::complex<double>(reference).imag(),
                                   tolerance),
             std::complex<double>(arg).imag(), result_listener);
}
} // namespace

namespace ae108 {
namespace cppslepc {
namespace {

template <class Policy> struct LinearEigenvalueProblemSolver_Test : Test {
  using matrix_type = cpppetsc::Matrix<Policy>;

  using solver_type = LinearEigenvalueProblemSolver<Policy>;
  solver_type solver;
};

using Policies =
    Types<cpppetsc::SequentialComputePolicy, cpppetsc::ParallelComputePolicy>;

TYPED_TEST_CASE(LinearEigenvalueProblemSolver_Test, Policies);

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           finds_the_two_eigenvalues_of_diagonal_matrix) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  EXPECT_THAT(this->solver.solve(A),
              ElementsAre(ComplexNear(2., 1e-7), ComplexNear(1., 1e-7)));
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test, solve_can_be_called_twice) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  this->solver.solve(A);

  const auto B = matrix_type::fromList({
      {3., 0.},
      {0., 4.},
  });

  EXPECT_THAT(this->solver.solve(B),
              ElementsAre(ComplexNear(4., 1e-7), ComplexNear(3., 1e-7)));
}

TYPED_TEST(
    LinearEigenvalueProblemSolver_Test,
    generalized_eigenvalue_problem_with_identity_matrix_returns_eigenvalues) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  const auto B = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  EXPECT_THAT(this->solver.solve(A, B),
              ElementsAre(ComplexNear(2., 1e-7), ComplexNear(1., 1e-7)));
}

TYPED_TEST(
    LinearEigenvalueProblemSolver_Test,
    generalized_eigenvalue_problem_with_scaled_identity_matrix_returns_scaled_eigenvalues) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  const auto B = matrix_type::fromList({
      {2., 0.},
      {0., 2.},
  });

  EXPECT_THAT(
      this->solver.solve(A, B),
      ElementsAre(ComplexNear(2. / 2., 1e-7), ComplexNear(1. / 2., 1e-7)));
}

} // namespace
} // namespace cppslepc
} // namespace ae108