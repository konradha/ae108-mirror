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
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/solve/TAOSolver.h"
#include "ae108/solve/test/Solver_Test.h"
#include <gmock/gmock.h>
#include <mpi.h>
#include <petscsys.h>

using ae108::cppptest::ScalarNear;
using testing::Eq;
using testing::SizeIs;
using testing::Types;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct TAOSolverTestConfig {
  template <class Assembler> using solver_type = TAOSolver<Assembler>;
  using policy_type = Policy;
};
} // namespace

namespace test {
namespace {

using Configurations =
    Types<TAOSolverTestConfig<cpppetsc::SequentialComputePolicy>,
          TAOSolverTestConfig<cpppetsc::ParallelComputePolicy>>;
INSTANTIATE_TYPED_TEST_CASE_P(TAO, Solver_Test, Configurations);
} // namespace
} // namespace test

namespace {

struct TAOSolver_Test
    : test::Solver_Test<
          TAOSolverTestConfig<cpppetsc::SequentialComputePolicy>> {};

TEST_F(TAOSolver_Test, local_functional_interface_works) {
  auto guess =
      cpppetsc::tag<cpppetsc::DistributedTag>(vector_type::fromList({2., 3.}));

  const auto factor = 2.;
  const auto solution = this->solver.computeSolution(
      {}, std::move(guess), factor * assembler_type::constantTime,
      [&](const cpppetsc::local<vector_type> &input, const double time,
          real_type *const output) {
        this->assembler.assembleEnergy(input, time / factor, output);
      },
      [&](const cpppetsc::local<vector_type> &input, const double time,
          cpppetsc::local<vector_type> *const output) {
        this->assembler.assembleForceVector(input, time / factor, output);
      },
      [&](const cpppetsc::local<vector_type> &input, const double time,
          matrix_type *const output) {
        this->assembler.assembleStiffnessMatrix(input, time / factor, output);
      });
  EXPECT_THAT(solution.unwrap(), SizeIs(2));
  EXPECT_THAT(solution(0), ScalarNear(1., 1e-7));
  EXPECT_THAT(solution(1), ScalarNear(2., 1e-7));
}

TEST_F(TAOSolver_Test, distributed_functional_interface_works) {
  auto guess =
      cpppetsc::tag<cpppetsc::DistributedTag>(vector_type::fromList({2., 3.}));

  const auto factor = 2.;
  const auto solution = this->solver.computeSolution(
      {}, std::move(guess), factor * assembler_type::constantTime,
      [&](const cpppetsc::distributed<vector_type> &input, const double time,
          real_type *const output) {
        auto displacements = vector_type::fromLocalMesh(this->mesh);
        this->mesh.copyToLocalVector(input, &displacements);
        this->assembler.assembleEnergy(displacements, time / factor, output);
        ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, output, 1, MPIU_REAL, MPIU_SUM,
                                  policy_type::communicator()),
                    Eq(0));
      },
      [&](const cpppetsc::distributed<vector_type> &input, const double time,
          cpppetsc::distributed<vector_type> *const output) {
        auto displacements = vector_type::fromLocalMesh(this->mesh);
        this->mesh.copyToLocalVector(input, &displacements);
        auto forces = vector_type::fromLocalMesh(this->mesh);
        this->assembler.assembleForceVector(displacements, time / factor,
                                            &forces);
        this->mesh.addToGlobalVector(forces, output);
      },
      [&](const cpppetsc::distributed<vector_type> &input, const double time,
          matrix_type *const output) {
        auto displacements = vector_type::fromLocalMesh(this->mesh);
        this->mesh.copyToLocalVector(input, &displacements);
        this->assembler.assembleStiffnessMatrix(displacements, time / factor,
                                                output);
        output->finalize();
      });
  EXPECT_THAT(solution.unwrap(), SizeIs(2));
  EXPECT_THAT(solution(0), ScalarNear(1., 1e-7));
  EXPECT_THAT(solution(1), ScalarNear(2., 1e-7));
}

} // namespace
} // namespace solve
} // namespace ae108

#endif