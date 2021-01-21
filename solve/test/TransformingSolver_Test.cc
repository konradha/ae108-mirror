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

#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/createRhsTransform.h"
#include "ae108/cpppetsc/createTransformInput.h"
#include "ae108/cpppetsc/createTransformOutput.h"
#include "ae108/solve/TransformingSolver.h"
#include "ae108/solve/test/Solver_Test.h"
#include <gmock/gmock.h>

using testing::DoubleNear;
using testing::SizeIs;
using testing::Types;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct TransformingSolverTestConfig {
  template <class Assembler> using solver_type = TransformingSolver<Assembler>;
  using policy_type = Policy;
};

template <class Configuration>
struct TransformingSolver_Test : test::Solver_Test<Configuration> {};

using Configurations =
    Types<TransformingSolverTestConfig<cpppetsc::SequentialComputePolicy>,
          TransformingSolverTestConfig<cpppetsc::ParallelComputePolicy>>;
TYPED_TEST_CASE(TransformingSolver_Test, Configurations);

TYPED_TEST(TransformingSolver_Test, no_bc_solve_works) {
  using assembler_type = typename TestFixture::assembler_type;
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;
  using vector_type = typename TestFixture::vector_type;

  auto guess = vector_type::fromGlobalMesh(this->mesh);
  guess.unwrap().replace().elements({0, 1}, {2., 3.});

  const auto transform = [&]() {
    auto A = matrix_type::fromMesh(this->mesh);
    A.replaceRowsByEye({0, 1});
    auto b = vector_type::fromGlobalMesh(this->mesh);
    return typename solver_type::BoundaryConditionContainer{std::move(A),
                                                            std::move(b)};
  }();

  const auto solution = this->solver.computeSolution(
      transform, std::move(guess), assembler_type::constantTime,
      &this->assembler);

  auto fullSolution = vector_type::fromDistributed(apply(transform, solution));
  ;
  EXPECT_THAT(fullSolution(0), DoubleNear(1., 1e-6));
  EXPECT_THAT(fullSolution(1), DoubleNear(2., 1e-6));
}

TYPED_TEST(TransformingSolver_Test, bc_solve_works) {
  using assembler_type = typename TestFixture::assembler_type;
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;
  using vector_type = typename TestFixture::vector_type;

  const auto matrix = matrix_type::fromMesh(this->mesh);

  const auto transform = [&]() {
    auto A = createRhsTransform(matrix, 1);
    A.assemblyView().replace().elements({0, 1}, {0}, {0., 1.});
    auto b = createTransformOutput(A);
    b.unwrap().replace().elements({0, 1}, {-6., 0.});
    return typename solver_type::BoundaryConditionContainer{std::move(A),
                                                            std::move(b)};
  }();

  auto guess = createTransformInput(transform.matrix);
  guess.unwrap().replace().elements({0}, {2.});

  const auto solution = this->solver.computeSolution(
      transform, std::move(guess), assembler_type::constantTime,
      &this->assembler);

  auto fullSolution = vector_type::fromDistributed(apply(transform, solution));

  ASSERT_THAT(fullSolution.unwrap(), SizeIs(2));
  EXPECT_THAT(fullSolution(0), DoubleNear(-6., 1e-6));
  EXPECT_THAT(fullSolution(1), DoubleNear(1., 1e-6));
}

} // namespace
} // namespace solve
} // namespace ae108
