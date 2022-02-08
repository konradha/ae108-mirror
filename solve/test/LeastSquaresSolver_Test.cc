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
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/solve/LeastSquaresSolver.h"
#include "ae108/solve/test/Solver_Test.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarNear;
using testing::SizeIs;
using testing::Types;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct LeastSquaresSolverTestConfig {
  template <class Assembler> using solver_type = LeastSquaresSolver<Assembler>;
  using policy_type = Policy;
};

template <class Configuration>
struct LeastSquaresSolver_Test : test::Solver_Test<Configuration> {};

using Configurations =
    Types<LeastSquaresSolverTestConfig<cpppetsc::SequentialComputePolicy>,
          LeastSquaresSolverTestConfig<cpppetsc::ParallelComputePolicy>>;
TYPED_TEST_CASE(LeastSquaresSolver_Test, Configurations);

TYPED_TEST(LeastSquaresSolver_Test, solve_yields_correct_result) {
  using assembler_type = typename TestFixture::assembler_type;
  using vector_type = typename TestFixture::vector_type;

  auto guess = vector_type::fromGlobalMesh(this->mesh);

  const auto solution = this->solver.computeSolution(
      std::move(guess), assembler_type::constantTime, &this->assembler);

  const auto fullSolution = vector_type::fromDistributed(solution);

  ASSERT_THAT(fullSolution.unwrap(), SizeIs(2));
  EXPECT_THAT(fullSolution(0), ScalarNear(1., 1e-6));
  EXPECT_THAT(fullSolution(1), ScalarNear(2., 1e-6));
}

} // namespace
} // namespace solve
} // namespace ae108

#endif