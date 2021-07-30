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

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/solve/DynamicSolver.h"
#include "ae108/solve/NonlinearSolver.h"
#include "ae108/solve/dynamics/DynamicState.h"
#include "ae108/solve/dynamics/NewmarkParameters.h"
#include "ae108/solve/test/Solver_Test.h"
#include <array>
#include <functional>
#include <gmock/gmock.h>

using ae108::cppptest::AlmostEqIfLocal;
using testing::DoubleEq;
using testing::DoubleNear;
using testing::Eq;
using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct Solver_Mock {
  using reference_solver_type = NonlinearSolver<test::Assembler_Mock<Policy>>;

  using vector_type = typename reference_solver_type::vector_type;
  using distributed_vector_type = cpppetsc::distributed<vector_type>;
  using matrix_type = typename reference_solver_type::matrix_type;
  using mesh_type = typename reference_solver_type::mesh_type;

  using BoundaryConditionContainer =
      typename reference_solver_type::BoundaryConditionContainer;

  using DistributedForceVectorAssembler =
      typename reference_solver_type::DistributedForceVectorAssembler;
  using DistributedStiffnessMatrixAssembler =
      typename reference_solver_type::DistributedStiffnessMatrixAssembler;

  std::function<distributed_vector_type(
      const BoundaryConditionContainer &, distributed_vector_type, double,
      DistributedForceVectorAssembler, DistributedStiffnessMatrixAssembler)>
      computeSolution;
};

constexpr double constantTimeStep = .1;

template <class Policy> struct DynamicSolver_Test : Test {
  using assembler_type = test::Assembler_Mock<Policy>;
  assembler_type assembler;

  using solver_type = Solver_Mock<Policy>;
  solver_type solver;

  const dynamics::NewmarkParameters newmark = {.25, .5};

  using dynamic_solver_type = DynamicSolver<assembler_type, solver_type>;
  dynamic_solver_type dynamicSolver =
      dynamic_solver_type{&mesh, &solver, newmark};

  typename assembler_type::mesh_type mesh =
      assembler_type::mesh_type::template fromConnectivity<
          std::array<std::array<int, 1>, 1>>(1, {{{{0}}}}, 1, 2);

  using vector_type = typename dynamic_solver_type::vector_type;

  using matrix_type = typename dynamic_solver_type::matrix_type;
  matrix_type massMatrix = matrix_type::fromMesh(mesh);
  matrix_type dampingMatrix = matrix_type::fromMesh(mesh);

  using state_type = typename dynamic_solver_type::state_type;
  state_type state = {vector_type::fromGlobalMesh(mesh),
                      vector_type::fromGlobalMesh(mesh),
                      vector_type::fromGlobalMesh(mesh)};

  explicit DynamicSolver_Test() {
    state.displacements.unwrap().replace().elements({0, 1}, {2., 3.});
    state.velocities.unwrap().replace().elements({0, 1}, {1., 4.});
    state.accelerations.unwrap().replace().elements({0, 1}, {-1., -2.});
    massMatrix.assemblyView().replace().elements({0, 1}, {0, 1},
                                                 {-1., 2., -3., 4.});
    dampingMatrix.assemblyView().replace().elements({0, 1}, {0, 1},
                                                    {1., -3., 2., -4.});
  }
};

using Policies =
    Types<cpppetsc::SequentialComputePolicy, cpppetsc::ParallelComputePolicy>;
TYPED_TEST_CASE(DynamicSolver_Test, Policies);

TYPED_TEST(DynamicSolver_Test, passes_on_displacements) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          typename solver_type::distributed_vector_type guess, double,
          typename solver_type::DistributedForceVectorAssembler,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        const auto full = TestFixture::vector_type::fromDistributed(guess);
        EXPECT_THAT(full.unwrap(), SizeIs(2));
        EXPECT_THAT(full(0), DoubleEq(2.));
        EXPECT_THAT(full(1), DoubleEq(3.));
        return guess;
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test, passes_on_time) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          typename solver_type::distributed_vector_type guess,
          const double time,
          typename solver_type::DistributedForceVectorAssembler,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        EXPECT_THAT(time, DoubleEq(TestFixture::assembler_type::constantTime));
        return guess;
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test, passes_on_bc) {
  using solver_type = typename TestFixture::solver_type;
  typename solver_type::BoundaryConditionContainer boundaryConditions;
  for (const auto &vertex : this->mesh.localVertices()) {
    boundaryConditions.push_back({vertex, 1, .7});
  }

  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &bc,
          typename solver_type::distributed_vector_type guess, double,
          typename solver_type::DistributedForceVectorAssembler,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        EXPECT_THAT(bc, SizeIs(boundaryConditions.size()));
        if (!bc.empty()) {
          EXPECT_THAT(bc.front().node.index(),
                      Eq(boundaryConditions.front().node.index()));
          EXPECT_THAT(bc.front().dof, Eq(boundaryConditions.front().dof));
          EXPECT_THAT(bc.front().value, Eq(boundaryConditions.front().value));
        }
        return guess;
      };

  this->dynamicSolver.computeSolution(
      boundaryConditions, std::move(this->state),
      TestFixture::assembler_type::constantTime, constantTimeStep,
      this->massMatrix, this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test, returns_correct_displacements) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          typename solver_type::distributed_vector_type, double,
          typename solver_type::DistributedForceVectorAssembler,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        auto result = TestFixture::vector_type::fromGlobalMesh(this->mesh);
        result.unwrap().replace().elements({0, 1}, {7., 77.});
        return result;
      };

  const auto solution = this->dynamicSolver.computeSolution(
      {}, std::move(this->state), TestFixture::assembler_type::constantTime,
      constantTimeStep, this->massMatrix, this->dampingMatrix,
      &this->assembler);

  const auto displacements =
      TestFixture::vector_type::fromDistributed(solution.displacements);

  ASSERT_THAT(displacements.unwrap(), SizeIs(2));
  EXPECT_THAT(displacements(0), DoubleEq(7.));
  EXPECT_THAT(displacements(1), DoubleEq(77.));
}

TYPED_TEST(DynamicSolver_Test, returns_correct_velocities) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          typename solver_type::distributed_vector_type, double,
          typename solver_type::DistributedForceVectorAssembler,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        auto result = TestFixture::vector_type::fromGlobalMesh(this->mesh);
        result.unwrap().replace().elements({0, 1}, {7., 77.});
        return result;
      };

  const auto solution = this->dynamicSolver.computeSolution(
      {}, std::move(this->state), TestFixture::assembler_type::constantTime,
      constantTimeStep, this->massMatrix, this->dampingMatrix,
      &this->assembler);

  const auto velocities =
      TestFixture::vector_type::fromDistributed(solution.velocities);

  ASSERT_THAT(velocities.unwrap(), SizeIs(2));
  EXPECT_THAT(velocities(0), DoubleNear(99., 1e-7));
  EXPECT_THAT(velocities(1), DoubleNear(1476., 1e-7));
}

TYPED_TEST(DynamicSolver_Test, returns_correct_accelerations) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          typename solver_type::distributed_vector_type, double,
          typename solver_type::DistributedForceVectorAssembler,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        auto result = TestFixture::vector_type::fromGlobalMesh(this->mesh);
        result.unwrap().replace().elements({0, 1}, {7., 77.});
        return result;
      };

  const auto solution = this->dynamicSolver.computeSolution(
      {}, std::move(this->state), TestFixture::assembler_type::constantTime,
      constantTimeStep, this->massMatrix, this->dampingMatrix,
      &this->assembler);

  const auto accelerations =
      TestFixture::vector_type::fromDistributed(solution.accelerations);

  ASSERT_THAT(accelerations.unwrap(), SizeIs(2));
  EXPECT_THAT(accelerations(0), DoubleNear(1961., 1e-7));
  EXPECT_THAT(accelerations(1), DoubleNear(29442., 1e-7));
}

TYPED_TEST(DynamicSolver_Test, computes_correct_effective_forces) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          const typename solver_type::distributed_vector_type guess,
          const double time,
          const typename solver_type::DistributedForceVectorAssembler assemble,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        auto forces = TestFixture::vector_type::fromGlobalMesh(this->mesh);
        assemble(guess, time, &forces);

        const auto full = TestFixture::vector_type::fromDistributed(forces);

        EXPECT_THAT(full.unwrap(), SizeIs(2));
        EXPECT_THAT(full(0), DoubleEq(-260.));
        EXPECT_THAT(full(1), DoubleEq(-483.));

        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test,
           computes_same_effective_forces_when_called_twice) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          const typename solver_type::distributed_vector_type guess,
          const double time,
          const typename solver_type::DistributedForceVectorAssembler assemble,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        auto reference = TestFixture::vector_type::fromGlobalMesh(this->mesh);
        assemble(guess, time, &reference);

        auto forces = TestFixture::vector_type::fromGlobalMesh(this->mesh);
        assemble(guess, time, &forces);

        forces.unwrap().timesAlphaPlusBetaX(1., -1., reference);
        EXPECT_THAT(forces.unwrap().norm(), DoubleEq(0));

        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test,
           computes_correct_effective_forces_with_mod_probe) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          typename solver_type::distributed_vector_type guess,
          const double time,
          const typename solver_type::DistributedForceVectorAssembler assemble,
          typename solver_type::DistributedStiffnessMatrixAssembler) {
        guess.unwrap().scale(2.);
        const auto &probe = guess;

        auto forces = TestFixture::vector_type::fromGlobalMesh(this->mesh);
        assemble(probe, time, &forces);

        const auto full = TestFixture::vector_type::fromDistributed(forces);

        EXPECT_THAT(full.unwrap(), SizeIs(2));
        EXPECT_THAT(full(0), DoubleEq(1253.));
        EXPECT_THAT(full(1), DoubleEq(1944.));

        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test, computes_correct_effective_stiffness) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          const typename solver_type::distributed_vector_type guess,
          const double time,
          typename solver_type::DistributedForceVectorAssembler,
          const typename solver_type::DistributedStiffnessMatrixAssembler
              assemble) {
        auto matrix = TestFixture::matrix_type::fromMesh(this->mesh);
        assemble(guess, time, &matrix);

        EXPECT_THAT(matrix, AlmostEqIfLocal(0, 0, -368.));
        EXPECT_THAT(matrix, AlmostEqIfLocal(0, 1, 739.));
        EXPECT_THAT(matrix, AlmostEqIfLocal(1, 0, -1161.));
        EXPECT_THAT(matrix, AlmostEqIfLocal(1, 1, 1547.));

        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test,
           computes_same_effective_stiffness_when_called_twice) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          const typename solver_type::distributed_vector_type guess,
          const double time,
          typename solver_type::DistributedForceVectorAssembler,
          const typename solver_type::DistributedStiffnessMatrixAssembler
              assemble) {
        auto reference = TestFixture::matrix_type::fromMesh(this->mesh);
        assemble(guess, time, &reference);

        auto matrix = TestFixture::matrix_type::fromMesh(this->mesh);
        assemble(guess, time, &matrix);

        matrix.addAlphaX(-1., reference);
        EXPECT_THAT(matrix.norm(), DoubleEq(0.));

        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test,
           computes_correct_effective_stiffness_with_mod_probe) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [&](const typename solver_type::BoundaryConditionContainer &,
          typename solver_type::distributed_vector_type guess,
          const double time,
          typename solver_type::DistributedForceVectorAssembler,
          const typename solver_type::DistributedStiffnessMatrixAssembler
              assemble) {
        auto matrix = TestFixture::matrix_type::fromMesh(this->mesh);

        guess.unwrap().scale(2.);
        const auto &probe = guess;

        assemble(probe, time, &matrix);

        EXPECT_THAT(matrix, AlmostEqIfLocal(0, 0, -332.));
        EXPECT_THAT(matrix, AlmostEqIfLocal(0, 1, 739.));
        EXPECT_THAT(matrix, AlmostEqIfLocal(1, 0, -1161.));
        EXPECT_THAT(matrix, AlmostEqIfLocal(1, 1, 1628.));

        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  this->dynamicSolver.computeSolution({}, std::move(this->state),
                                      TestFixture::assembler_type::constantTime,
                                      constantTimeStep, this->massMatrix,
                                      this->dampingMatrix, &this->assembler);
}

TYPED_TEST(DynamicSolver_Test, local_functional_interface_is_callable) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [this](const typename solver_type::BoundaryConditionContainer &,
             const typename solver_type::distributed_vector_type, const double,
             typename solver_type::DistributedForceVectorAssembler,
             const typename solver_type::DistributedStiffnessMatrixAssembler) {
        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  const auto zeroMatrix = TestFixture::matrix_type::fromMesh(this->mesh);

  this->dynamicSolver.computeSolution(
      {}, std::move(this->state), TestFixture::assembler_type::constantTime,
      constantTimeStep, zeroMatrix, zeroMatrix,
      [](const cpppetsc::local<typename TestFixture::vector_type> &, double,
         cpppetsc::local<typename TestFixture::vector_type> *output) {
        output->unwrap().setZero();
      },
      [](const cpppetsc::local<typename TestFixture::vector_type> &, double,
         typename TestFixture::matrix_type *output) { output->setZero(); });
}

TYPED_TEST(DynamicSolver_Test, distributed_functional_interface_is_callable) {
  using solver_type = typename TestFixture::solver_type;
  this->solver.computeSolution =
      [this](const typename solver_type::BoundaryConditionContainer &,
             const typename solver_type::distributed_vector_type, const double,
             typename solver_type::DistributedForceVectorAssembler,
             const typename solver_type::DistributedStiffnessMatrixAssembler) {
        return TestFixture::vector_type::fromGlobalMesh(this->mesh);
      };

  const auto zeroMatrix = TestFixture::matrix_type::fromMesh(this->mesh);

  this->dynamicSolver.computeSolution(
      {}, std::move(this->state), TestFixture::assembler_type::constantTime,
      constantTimeStep, zeroMatrix, zeroMatrix,
      [](const cpppetsc::distributed<typename TestFixture::vector_type> &,
         double,
         cpppetsc::distributed<typename TestFixture::vector_type> *output) {
        output->unwrap().setZero();
      },
      [](const cpppetsc::distributed<typename TestFixture::vector_type> &,
         double,
         typename TestFixture::matrix_type *output) { output->setZero(); });
}
} // namespace
} // namespace solve
} // namespace ae108
