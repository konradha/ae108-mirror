// © 2022 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/cppptest/Matchers.h"
#include "ae108/solve/GeneralizedTAOSolver.h"
#include "ae108/solve/test/Solver_Test.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::Types;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct Configuration {
  using policy_type = Policy;
  template <class T> using solver_type = GeneralizedTAOSolver<T>;
};

template <class Configuration>
struct GeneralizedTAOSolver_Test : test::Solver_Test<Configuration> {};

using Configurations = Types<Configuration<cpppetsc::SequentialComputePolicy>,
                             Configuration<cpppetsc::ParallelComputePolicy>>;
TYPED_TEST_CASE(GeneralizedTAOSolver_Test, Configurations);

TYPED_TEST(GeneralizedTAOSolver_Test, solves_unconstrained_problem) {
  using assembler_type = typename TestFixture::assembler_type;
  using vector_type = typename TestFixture::vector_type;

  const auto solution = this->solver.computeSolution(
      {}, vector_type::fromGlobalMesh(this->mesh), assembler_type::constantTime,
      &this->assembler);

  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 1.));
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(1, 2.));
}

TYPED_TEST(GeneralizedTAOSolver_Test, solves_partially_constrained_problem) {
  using assembler_type = typename TestFixture::assembler_type;
  using matrix_type = typename TestFixture::matrix_type;
  using vector_type = typename TestFixture::vector_type;

  const auto solution = this->solver.computeSolution(
      {
          1,
          [](const cpppetsc::distributed<vector_type> &in,
             cpppetsc::distributed<vector_type> *const out) {
            const auto full = vector_type::fromDistributed(in);
            out->unwrap().setZero();
            out->unwrap().replace().element(0, full(0) - full(1));
          },
          [](const cpppetsc::distributed<vector_type> &,
             matrix_type *const out) {
            out->setZero();
            out->preallocatedAssemblyView(2)
                .replace()
                .element(0, 0, 1.)
                .element(0, 1, -1.);
          },
      },
      vector_type::fromGlobalMesh(this->mesh), assembler_type::constantTime,
      &this->assembler);

  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 3. / 2.));
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(1, 3. / 2.));
}

TYPED_TEST(GeneralizedTAOSolver_Test, solves_fully_constrained_problem) {
  using assembler_type = typename TestFixture::assembler_type;
  using matrix_type = typename TestFixture::matrix_type;
  using vector_type = typename TestFixture::vector_type;

  const auto solution = this->solver.computeSolution(
      {
          2,
          [](const cpppetsc::distributed<vector_type> &in,
             cpppetsc::distributed<vector_type> *const out) {
            const auto full = vector_type::fromDistributed(in);
            out->unwrap().setZero();
            out->unwrap()
                .replace()
                .element(0, full(0) - 3.)
                .element(1, full(1) - 4.);
          },
          [](const cpppetsc::distributed<vector_type> &,
             matrix_type *const out) {
            out->setZero();
            out->preallocatedAssemblyView(1)
                .replace()
                .element(0, 0, 1.)
                .element(1, 1, 1.);
          },
      },
      vector_type::fromGlobalMesh(this->mesh), assembler_type::constantTime,
      &this->assembler);

  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 3.));
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(1, 4.));
}

TYPED_TEST(GeneralizedTAOSolver_Test,
           assembler_can_be_replaced_by_local_assembly_functions) {
  using solver_type = typename TestFixture::solver_type;
  using matrix_type = typename TestFixture::matrix_type;
  using vector_type = typename TestFixture::vector_type;

  const auto solution = this->solver.computeSolution(
      {}, vector_type::fromGlobalMesh(this->mesh),
      TestFixture::assembler_type::constantTime,
      [assembler = &this->assembler](
          const cpppetsc::local<vector_type> &in, const double time,
          typename solver_type::real_type *const out) {
        assembler->assembleEnergy(in, time, out);
      },
      [assembler = &this->assembler](const cpppetsc::local<vector_type> &in,
                                     const double time,
                                     cpppetsc::local<vector_type> *const out) {
        assembler->assembleForceVector(in, time, out);
      },
      [assembler = &this->assembler](const cpppetsc::local<vector_type> &in,
                                     const double time,
                                     matrix_type *const out) {
        assembler->assembleStiffnessMatrix(in, time, out);
      });

  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 1.));
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(1, 2.));
}

TYPED_TEST(GeneralizedTAOSolver_Test,
           assembler_can_be_replaced_by_distributed_assembly_functions) {
  using policy_type = typename TestFixture::policy_type;
  using vector_type = typename TestFixture::vector_type;
  using matrix_type = typename TestFixture::matrix_type;
  using real_type = typename TestFixture::solver_type::real_type;

  const auto solution = this->solver.computeSolution(
      {}, vector_type::fromGlobalMesh(this->mesh),
      TestFixture::assembler_type::constantTime,
      [assembler = &this->assembler,
       mesh = &this->mesh](const cpppetsc::distributed<vector_type> &in,
                           const double time, real_type *const out) {
        auto local = vector_type::fromLocalMesh(*mesh);
        mesh->copyToLocalVector(in, &local);
        assembler->assembleEnergy(local, time, out);
        policy_type::handleError(MPI_Allreduce(MPI_IN_PLACE, out, 1, MPIU_REAL,
                                               MPIU_SUM,
                                               policy_type::communicator()));
      },
      [assembler = &this->assembler, mesh = &this->mesh](
          const cpppetsc::distributed<vector_type> &in, const double time,
          cpppetsc::distributed<vector_type> *const out) {
        auto localIn = vector_type::fromLocalMesh(*mesh);
        auto localOut = vector_type::fromLocalMesh(*mesh);
        mesh->copyToLocalVector(in, &localIn);
        assembler->assembleForceVector(localIn, time, &localOut);
        mesh->addToGlobalVector(localOut, out);
      },
      [assembler = &this->assembler,
       mesh = &this->mesh](const cpppetsc::distributed<vector_type> &in,
                           const double time, matrix_type *const out) {
        auto local = vector_type::fromLocalMesh(*mesh);
        mesh->copyToLocalVector(in, &local);
        assembler->assembleStiffnessMatrix(local, time, out);
        out->finalize();
      });

  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(0, 1.));
  EXPECT_THAT(solution.unwrap(), ScalarEqIfLocal(1, 2.));
}
} // namespace
} // namespace solve
} // namespace ae108

#endif