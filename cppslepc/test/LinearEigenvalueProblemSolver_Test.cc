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

#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/scalarProduct.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/cppslepc/InvalidEigenvalueIndexException.h"
#include "ae108/cppslepc/InvalidProblemTypeException.h"
#include "ae108/cppslepc/LinearEigenvalueProblemSolver.h"
#include "ae108/cppslepc/NoOperatorsSetException.h"
#include <gmock/gmock.h>

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
  using solver_type = LinearEigenvalueProblemSolver<Policy>;
  using vector_type = typename solver_type::vector_type;
  using matrix_type = typename solver_type::matrix_type;
  using eigenpair_type = EigenPair<Policy>;
};

using Policies =
    Types<cpppetsc::SequentialComputePolicy, cpppetsc::ParallelComputePolicy>;
TYPED_TEST_CASE(LinearEigenvalueProblemSolver_Test, Policies);

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           hermitian_evp_eigen_pair_is_correct) {
  using solver_type = typename TestFixture::solver_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;
  using eigenpair_type = typename TestFixture::eigenpair_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  solver.setOperators(&A, solver_type::Type::hermitian);
  solver.solve();

  ASSERT_THAT(solver.numberOfEigenpairs(), Eq(2));

#ifdef AE108_PETSC_COMPLEX
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2)};
#else
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2), vector_type(2)};
#endif

  solver.getEigenpair(0, &eigenpair);
  EXPECT_THAT(eigenpair.value, ComplexNear(2, 1e-15));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(
      scalarProduct(eigenpair.vector,
                    vector_type(vector_type::value_type::fromList({1., 0.}))),
      cppptest::ScalarNear(0., 1e-15));
#else
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_real,
                    vector_type(vector_type::value_type::fromList({1., 0.}))),
      cppptest::ScalarNear(0., 1e-15));
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_imag,
                    vector_type(vector_type::value_type::fromList({1., 0.}))),
      cppptest::ScalarNear(0., 1e-15));
#endif

  solver.getEigenpair(1, &eigenpair);

  EXPECT_THAT(eigenpair.value, ComplexNear(1, 1e-15));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(
      scalarProduct(eigenpair.vector,
                    vector_type(vector_type::value_type::fromList({0., 1.}))),
      cppptest::ScalarNear(0., 1e-15));
#else
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_real,
                    vector_type(vector_type::value_type::fromList({0., 1.}))),
      cppptest::ScalarNear(0., 1e-15));
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_imag,
                    vector_type(vector_type::value_type::fromList({0., 1.}))),
      cppptest::ScalarNear(0., 1e-15));
#endif
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           generalized_evp_eigen_pair_is_correct) {
  using solver_type = typename TestFixture::solver_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;
  using eigenpair_type = typename TestFixture::eigenpair_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  const auto B = matrix_type::fromList({
      {3., 0.},
      {0., 3.},
  });

  solver.setOperators(&A, &B);
  solver.solve();

  ASSERT_THAT(solver.numberOfEigenpairs(), Eq(2));

#ifdef AE108_PETSC_COMPLEX
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2)};
#else
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2), vector_type(2)};
#endif

  solver.getEigenpair(0, &eigenpair);
  EXPECT_THAT(eigenpair.value, ComplexNear(2. / 3., 1e-15));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(
      scalarProduct(eigenpair.vector,
                    vector_type(vector_type::value_type::fromList({1., 0.}))),
      cppptest::ScalarNear(0., 1e-15));
#else
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_real,
                    vector_type(vector_type::value_type::fromList({1., 0.}))),
      cppptest::ScalarNear(0., 1e-15));
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_imag,
                    vector_type(vector_type::value_type::fromList({1., 0.}))),
      cppptest::ScalarNear(0., 1e-15));
#endif

  solver.getEigenpair(1, &eigenpair);

  EXPECT_THAT(eigenpair.value, ComplexNear(1. / 3., 1e-15));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(
      scalarProduct(eigenpair.vector,
                    vector_type(vector_type::value_type::fromList({0., 1.}))),
      cppptest::ScalarNear(0., 1e-15));
#else
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_real,
                    vector_type(vector_type::value_type::fromList({0., 1.}))),
      cppptest::ScalarNear(0., 1e-15));
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_imag,
                    vector_type(vector_type::value_type::fromList({0., 1.}))),
      cppptest::ScalarNear(0., 1e-15));
#endif
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           nonhermitian_evp_eigen_pair_is_correct) {
  using solver_type = typename TestFixture::solver_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;
  using eigenpair_type = typename TestFixture::eigenpair_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 4.},
      {1., 1.},
  });

  solver.setOperators(&A, solver_type::Type::nonhermitian);
  solver.solve();

  ASSERT_THAT(solver.numberOfEigenpairs(), Eq(2));

#ifdef AE108_PETSC_COMPLEX
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2)};
#else
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2), vector_type(2)};
#endif

  solver.getEigenpair(0, &eigenpair);
  EXPECT_THAT(eigenpair.value, ComplexNear(3, 1e-15));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(
      scalarProduct(eigenpair.vector,
                    vector_type(vector_type::value_type::fromList({1., -2.}))),
      cppptest::ScalarNear(0., 1e-15));
#else
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_real,
                    vector_type(vector_type::value_type::fromList({1., -2.}))),
      cppptest::ScalarNear(0., 1e-15));
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_imag,
                    vector_type(vector_type::value_type::fromList({1., -2.}))),
      cppptest::ScalarNear(0., 1e-15));
#endif

  solver.getEigenpair(1, &eigenpair);

  EXPECT_THAT(eigenpair.value, ComplexNear(-1, 1e-15));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(
      scalarProduct(eigenpair.vector,
                    vector_type(vector_type::value_type::fromList({1., 2.}))),
      cppptest::ScalarNear(0., 1e-15));
#else
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_real,
                    vector_type(vector_type::value_type::fromList({1., 2.}))),
      cppptest::ScalarNear(0., 1e-15));
  EXPECT_THAT(
      scalarProduct(eigenpair.vector_imag,
                    vector_type(vector_type::value_type::fromList({1., 2.}))),
      cppptest::ScalarNear(0., 1e-15));
#endif
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           throws_exception_for_nonexisting_eigenvalue_index) {
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  solver.setOperators(&A);

  ASSERT_THAT(solver.numberOfEigenpairs(), Eq(0));
  EXPECT_THROW(solver.getEigenvalue(-1), InvalidEigenvalueIndexException);
  EXPECT_THROW(solver.getEigenvalue(0), InvalidEigenvalueIndexException);
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           throws_exception_for_nonexisting_eigenpair_index) {
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;
  using vector_type = typename TestFixture::vector_type;
  using eigenpair_type = typename TestFixture::eigenpair_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  solver.setOperators(&A);

#ifdef AE108_PETSC_COMPLEX
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2)};
#else
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2), vector_type(2)};
#endif

  ASSERT_THAT(solver.numberOfEigenpairs(), Eq(0));
  EXPECT_THROW(solver.getEigenpair(-1, &eigenpair),
               InvalidEigenvalueIndexException);
  EXPECT_THROW(solver.getEigenpair(0, &eigenpair),
               InvalidEigenvalueIndexException);
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           throws_exception_if_operators_not_set) {
  using solver_type = typename TestFixture::solver_type;

  auto solver = solver_type{};

  EXPECT_THROW(solver.solve(), NoOperatorsSetException);
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           throws_exception_if_generalized_type_for_one_operator) {
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  EXPECT_THROW(
      solver.setOperators(&A, solver_type::Type::generalized_hermitian),
      InvalidProblemTypeException);
  EXPECT_THROW(
      solver.setOperators(&A, solver_type::Type::generalized_nonhermitian),
      InvalidProblemTypeException);
  EXPECT_THROW(
      solver.setOperators(&A, solver_type::Type::generalized_nonhermitian_spd),
      InvalidProblemTypeException);
  EXPECT_THROW(
      solver.setOperators(&A, solver_type::Type::generalized_indefinite),
      InvalidProblemTypeException);
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           throws_exception_if_nongeneralized_type_for_two_operators) {
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });

  EXPECT_THROW(solver.setOperators(&A, &A, solver_type::Type::hermitian),
               InvalidProblemTypeException);
  EXPECT_THROW(solver.setOperators(&A, &A, solver_type::Type::nonhermitian),
               InvalidProblemTypeException);
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           solver_supports_switching_problem_type_to_hermitian_type) {
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });
  const auto B = matrix_type::fromList({
      {2., 0.},
      {0., 2.},
  });

  solver.setOperators(&A, &B, solver_type::Type::generalized_hermitian);
  solver.setOperators(&A, solver_type::Type::hermitian);

  solver.solve();

  EXPECT_THAT(solver.getEigenvalue(0), cppptest::ScalarNear(1., 1e-15));
  EXPECT_THAT(solver.getEigenvalue(1), cppptest::ScalarNear(1., 1e-15));
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           solver_supports_switching_problem_type_to_unspecified_type) {
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;

  auto solver = solver_type{};

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 1.},
  });
  const auto B = matrix_type::fromList({
      {2., 0.},
      {0., 2.},
  });

  solver.setOperators(&A, &B, solver_type::Type::generalized_hermitian);
  solver.setOperators(&A);

  solver.solve();

  EXPECT_THAT(solver.getEigenvalue(0), cppptest::ScalarNear(1., 1e-15));
  EXPECT_THAT(solver.getEigenvalue(1), cppptest::ScalarNear(1., 1e-15));
}

} // namespace
} // namespace cppslepc
} // namespace ae108
