// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/asThAT.h"
#include "ae108/cpppetsc/computeElementsOfMatrix.h"
#include "ae108/cppslepc/computeSmallestEigenvalues.h"
#include <gmock/gmock.h>

using testing::ElementsAre;
using testing::Eq;
using testing::Ge;
using testing::IsEmpty;
using testing::SizeIs;
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

template <class Policy> struct computeSmallestEigenvalues_Test : Test {
  using matrix_type = cpppetsc::Matrix<Policy>;
};

using Policies =
    Types<cpppetsc::SequentialComputePolicy, cpppetsc::ParallelComputePolicy>;

TYPED_TEST_CASE(computeSmallestEigenvalues_Test, Policies);

TYPED_TEST(computeSmallestEigenvalues_Test,
           finds_the_two_eigenvalues_of_diagonal_matrix) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  EXPECT_THAT(computeSmallestEigenvalues(A, 2),
              ElementsAre(ComplexNear(1., 1e-7), ComplexNear(2., 1e-7)));
}

TYPED_TEST(computeSmallestEigenvalues_Test,
           finds_at_least_smallest_eigenvalue_of_diagonal_matrix) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  const auto result = computeSmallestEigenvalues(A, 1);
  ASSERT_THAT(result, SizeIs(Ge(1)));
  EXPECT_THAT(result.front(), ComplexNear(1., 1e-7));
}

TYPED_TEST(computeSmallestEigenvalues_Test,
           finding_zero_eigenvalues_returns_empty_list) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  EXPECT_THAT(computeSmallestEigenvalues(A, 0), IsEmpty());
}

TYPED_TEST(computeSmallestEigenvalues_Test,
           requesting_too_many_eigenvalues_returns_all_eigenvalues) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  EXPECT_THAT(computeSmallestEigenvalues(A, 2),
              ElementsAre(ComplexNear(1., 1e-7), ComplexNear(2., 1e-7)));
}

TYPED_TEST(computeSmallestEigenvalues_Test,
           finding_zero_generalized_eigenvalues_returns_empty_list) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  const auto B = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  EXPECT_THAT(computeSmallestEigenvalues(A, B, 0), IsEmpty());
}

TYPED_TEST(computeSmallestEigenvalues_Test,
           finds_the_two_eigenvalues_of_transformed_diagonal_matrix) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  const auto T = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  EXPECT_THAT(computeSmallestEigenvalues(cpppetsc::asThAT(&A, &T), 2),
              ElementsAre(ComplexNear(1., 1e-7), ComplexNear(2., 1e-7)));
}

TYPED_TEST(
    computeSmallestEigenvalues_Test,
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

  EXPECT_THAT(computeSmallestEigenvalues(A, B, 2),
              ElementsAre(ComplexNear(1., 1e-7), ComplexNear(2., 1e-7)));
}

TYPED_TEST(
    computeSmallestEigenvalues_Test,
    transformed_generalized_eigenvalue_problem_with_identity_matrix_returns_eigenvalues) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });
  const auto B = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });
  const auto T = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  EXPECT_THAT(computeSmallestEigenvalues(
                  cpppetsc::asThAT(&A, &T),
                  cpppetsc::computeElementsOfMatrix(cpppetsc::asThAT(&B, &T)),
                  2),
              ElementsAre(ComplexNear(1., 1e-7), ComplexNear(2., 1e-7)));
}

TYPED_TEST(
    computeSmallestEigenvalues_Test,
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
      computeSmallestEigenvalues(A, B, 2),
      ElementsAre(ComplexNear(1. / 2., 1e-7), ComplexNear(2. / 2., 1e-7)));
}

} // namespace
} // namespace cppslepc
} // namespace ae108