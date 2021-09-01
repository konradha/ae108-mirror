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

#pragma once

#include "ae108/assembly/AssemblerTypeTraits.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cppptest/Matchers.h"
#include <array>
#include <gmock/gmock.h>

namespace ae108 {
namespace solve {
namespace test {
namespace {

template <class Policy> struct Assembler_Mock {
  static constexpr double constantTime = .77;

  using policy_type = Policy;
  using mesh_type = cpppetsc::Mesh<policy_type>;
  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;
  using real_type = typename mesh_type::real_type;

  /**
   * @brief Computes |x - 1|^2 * |y - 2|^2 + 1
   * @remark Expects to be called with the time constant.
   * @remark Expects two degrees of freedom.
   */
  void assembleEnergy(const cpppetsc::local<vector_type> &displacements,
                      const double time, real_type *const energy) const {
    assert(energy);

    EXPECT_THAT(time, ::testing::DoubleEq(constantTime));

    if (displacements.unwrap().size() > 0) {
      const auto x = displacements(0);
      const auto y = displacements(1);

      *energy = std::norm(x - 1) + std::norm(y - 2) + 1.;
    }
  }

  /**
   * @brief Computes the derivative of the energy.
   * @remark Expects to be called with the time constant.
   * @remark Expects two degrees of freedom.
   */
  void assembleForceVector(const cpppetsc::local<vector_type> &displacements,
                           const double time,
                           cpppetsc::local<vector_type> *const force) const {
    assert(force);

    EXPECT_THAT(time, ::testing::DoubleEq(constantTime));

    const auto replacer = force->unwrap().replace();
    if (displacements.unwrap().size() > 0) {
      const auto x = displacements(0);
      const auto y = displacements(1);

      replacer(0) = 2. * std::real(x - 1.);
      replacer(1) = 2. * std::real(y - 2.);
    }
  }

  /**
   * @brief Computes the derivative of the forces.
   * @remark Expects to be called with the time constant.
   * @remark Expects two degrees of freedom, both on the same rank.
   */
  void
  assembleStiffnessMatrix(const cpppetsc::local<vector_type> &displacements,
                          const double time, matrix_type *const matrix) const {
    assert(matrix);

    EXPECT_THAT(time, ::testing::DoubleEq(constantTime));

    const auto replacer = matrix->assemblyView().replace();
    if (displacements.unwrap().size() > 0) {
      const auto x = displacements(0);
      const auto y = displacements(1);

      replacer(0, 0) = 2.;
      replacer(0, 1) = 0.;
      replacer(1, 0) = 0.;
      replacer(1, 1) = 2.;
    }
  }
};
} // namespace
} // namespace test
} // namespace solve
} // namespace ae108

namespace ae108 {
namespace assembly {

template <class Policy>
struct PolicyTypeTrait<solve::test::Assembler_Mock<Policy>> {
  using type = Policy;
};

template <class Policy>
struct MeshTypeTrait<solve::test::Assembler_Mock<Policy>> {
  using type = typename solve::test::Assembler_Mock<Policy>::mesh_type;
};
} // namespace assembly
} // namespace ae108

namespace ae108 {
namespace solve {
namespace test {
namespace {

template <class TestConfiguration> struct Solver_Test : ::testing::Test {
  using policy_type = typename TestConfiguration::policy_type;
  using assembler_type = Assembler_Mock<policy_type>;
  using mesh_type = typename assembler_type::mesh_type;
  using value_type = typename assembler_type::value_type;
  using real_type = typename assembler_type::real_type;
  using vector_type = typename assembler_type::vector_type;
  using matrix_type = typename assembler_type::matrix_type;
  using solver_type =
      typename TestConfiguration::template solver_type<assembler_type>;

  mesh_type mesh =
      mesh_type::template fromConnectivity<std::array<std::array<int, 2>, 2>>(
          1, {{{0, 1}, {0, 1}}}, 2, 1);

  solver_type solver = solver_type{&mesh};
  assembler_type assembler;
};

TYPED_TEST_CASE_P(Solver_Test);

TYPED_TEST_P(Solver_Test, no_bc_solve_works) {
  auto guess = TestFixture::vector_type::fromGlobalMesh(this->mesh);
  guess.unwrap().replace().elements({0, 1}, {2., 3.});

  const auto solution = this->solver.computeSolution(
      {}, std::move(guess), TestFixture::assembler_type::constantTime,
      &this->assembler);

  const auto fullSolution = TestFixture::vector_type::fromDistributed(solution);
  EXPECT_THAT(fullSolution.unwrap(), ::testing::SizeIs(2));
  EXPECT_THAT(fullSolution(0), ::ae108::cppptest::ValueNear(1., 1e-6));
  EXPECT_THAT(fullSolution(1), ::ae108::cppptest::ValueNear(2., 1e-6));
}

TYPED_TEST_P(Solver_Test, bc_solve_works) {
  auto guess = TestFixture::vector_type::fromGlobalMesh(this->mesh);
  guess.unwrap().replace().elements({0, 1}, {2., 3.});

  auto bc = typename TestFixture::solver_type::BoundaryConditionContainer{};
  for (const auto &vertex : this->mesh.localVertices()) {
    const auto range = vertex.localDofLineRange();
    if (range.first <= 0 && 0 < range.second) {
      bc.push_back({vertex, 0, -6.});
    }
  }

  const auto solution = this->solver.computeSolution(
      bc, std::move(guess), TestFixture::assembler_type::constantTime,
      &this->assembler);
  const auto fullSolution = TestFixture::vector_type::fromDistributed(solution);
  EXPECT_THAT(fullSolution.unwrap(), ::testing::SizeIs(2));
  EXPECT_THAT(fullSolution(0), ::ae108::cppptest::ValueNear(-6., 1e-6));
  EXPECT_THAT(fullSolution(1), ::ae108::cppptest::ValueNear(2., 1e-6));
}

TYPED_TEST_P(Solver_Test, full_bc_solve_works) {
  auto guess = TestFixture::vector_type::fromGlobalMesh(this->mesh);
  guess.unwrap().replace().elements({0, 1}, {2., 3.});

  auto bc = typename TestFixture::solver_type::BoundaryConditionContainer{};
  for (const auto &vertex : this->mesh.localVertices()) {
    const auto range = vertex.localDofLineRange();
    if (range.first <= 0 && 0 < range.second) {
      bc.push_back({vertex, 0, 7.});
    }
    if (range.first <= 1 && 1 < range.second) {
      bc.push_back({vertex, 0, 77.});
    }
  }

  const auto solution = this->solver.computeSolution(
      bc, std::move(guess), TestFixture::assembler_type::constantTime,
      &this->assembler);
  const auto fullSolution = TestFixture::vector_type::fromDistributed(solution);
  EXPECT_THAT(fullSolution.unwrap(), ::testing::SizeIs(2));
  EXPECT_THAT(fullSolution(0), ::ae108::cppptest::ValueNear(7., 1e-6));
  EXPECT_THAT(fullSolution(1), ::ae108::cppptest::ValueNear(77., 1e-6));
}

REGISTER_TYPED_TEST_CASE_P(Solver_Test, bc_solve_works, no_bc_solve_works,
                           full_bc_solve_works);

} // namespace
} // namespace test
} // namespace solve
} // namespace ae108
