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