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
#include "ae108/cppptest/Matchers.h"
#include "ae108/cppslepc/LinearEigenvalueProblemSolver.h"
#include <gmock/gmock.h>

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
           hermetian_evp_eigen_pair_is_correct) {
  using solver_type = typename TestFixture::solver_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;
  using eigenpair_type = typename TestFixture::eigenpair_type;

  auto solver = solver_type{};
  EPSSetProblemType(solver.data(), EPS_HEP);

  const auto A = matrix_type::fromList({
      {1., 0.},
      {0., 2.},
  });

  solver.setOperators(&A);

  solver.solve();

#ifdef AE108_PETSC_COMPLEX
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2)};
#else
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2), vector_type(2)};
#endif

  solver.getEigenpair(0, &eigenpair);
  EXPECT_THAT(eigenpair.value, ComplexNear(2, 1e-7));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(eigenpair.vector(0) * 1. + eigenpair.vector(1) * 0,
              cppptest::ScalarNear(0., 1e-7));
#else
  EXPECT_THAT(eigenpair.vector_real(0) * 1. + eigenpair.vector_real(1) * 0,
              cppptest::ScalarNear(0., 1e-7));
  EXPECT_THAT(eigenpair.vector_imag(0) * 0. + eigenpair.vector_imag(1) * 0,
              cppptest::ScalarNear(0., 1e-7));
#endif

  solver.getEigenpair(1, &eigenpair);

  EXPECT_THAT(eigenpair.value, ComplexNear(1, 1e-7));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT(eigenpair.vector(0) * 0. + eigenpair.vector(1) * 1,
              cppptest::ScalarNear(0., 1e-7));
#else
  EXPECT_THAT(eigenpair.vector_real(0) * 0. + eigenpair.vector_real(1) * 1,
              cppptest::ScalarNear(0., 1e-7));
  EXPECT_THAT(eigenpair.vector_imag(0) * 0. + eigenpair.vector_imag(1) * 0,
              cppptest::ScalarNear(0., 1e-7));
#endif
}

TYPED_TEST(LinearEigenvalueProblemSolver_Test,
           non_hermetian_evp_eigen_pair_is_correct) {
  using solver_type = typename TestFixture::solver_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;
  using eigenpair_type = typename TestFixture::eigenpair_type;

  auto solver = solver_type{};
  EPSSetProblemType(solver.data(), EPS_NHEP);

  const auto A = matrix_type::fromList({
      {1., -1.},
      {1., 1.},
  });

  solver.setOperators(&A);

  solver.solve();

#ifdef AE108_PETSC_COMPLEX
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2)};
#else
  auto eigenpair = eigenpair_type{{0, 0}, vector_type(2), vector_type(2)};
#endif

  solver.getEigenpair(0, &eigenpair);
  EXPECT_THAT(eigenpair.value, ComplexNear(std::complex<double>{1, 1}, 1e-7));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT((eigenpair.vector(0) * std::complex<double>{1, 0} +
               eigenpair.vector(1) * std::complex<double>{0, 1}),
              ComplexNear(0., 1e-7));
#else
  EXPECT_THAT(eigenpair.vector_real(0) * 1. + eigenpair.vector_real(1) * 0,
              cppptest::ScalarNear(0., 1e-7));
  EXPECT_THAT(eigenpair.vector_imag(0) * 0. + eigenpair.vector_imag(1) * 1.,
              cppptest::ScalarNear(0., 1e-7));
#endif

  solver.getEigenpair(1, &eigenpair);

  EXPECT_THAT(eigenpair.value, ComplexNear(std::complex<double>{1, -1}, 1e-7));
#ifdef AE108_PETSC_COMPLEX
  EXPECT_THAT((eigenpair.vector(0) * std::complex<double>{1, 0} +
               eigenpair.vector(1) * std::complex<double>{0, -1}),
              ComplexNear(0., 1e-7));
#else
  EXPECT_THAT(eigenpair.vector_real(0) * 1. + eigenpair.vector_real(1) * 0.,
              cppptest::ScalarNear(0., 1e-7));
  EXPECT_THAT(eigenpair.vector_imag(0) * 0. + eigenpair.vector_imag(1) * -1.,
              cppptest::ScalarNear(0., 1e-7));
#endif
}

} // namespace
} // namespace cppslepc
} // namespace ae108
