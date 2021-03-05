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

#include "ae108/assembly/Assembler.h"
#include "ae108/assembly/AssemblerTypeTraits.h"
#include "ae108/assembly/FeaturePlugin.h"
#include "ae108/assembly/FeaturePlugins.h"
#include "ae108/assembly/plugins/AssembleConsistentMassMatrixPlugin.h"
#include "ae108/assembly/plugins/AssembleEnergyPlugin.h"
#include "ae108/assembly/plugins/AssembleForceIfPlugin.h"
#include "ae108/assembly/plugins/AssembleForceVectorPlugin.h"
#include "ae108/assembly/plugins/AssembleLumpedMassMatrixPlugin.h"
#include "ae108/assembly/plugins/AssembleStiffnessMatrixPlugin.h"
#include "ae108/assembly/plugins/UpdateInternalVariablesPlugin.h"
#include "ae108/assembly/test/Element.h"
#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>
#include <mpi.h>
#include <petscsys.h>
#include <string>
#include <type_traits>
#include <vector>

using ae108::cppptest::AlmostEqIfLocal;
using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::NiceMock;
using testing::Pair;
using testing::PrintToString;
using testing::Return;
using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace assembly {
namespace {

DEFINE_CONST_ASSEMBLER_PLUGIN(
    SumUpElementIndicesPlugin, sumUpElementIndices,
    (typename SizeTypeTrait<Assembler>::type *const output)) {
  for (const auto &annotatedElement : this->assembler().meshElements()) {
    *output += annotatedElement.meshView().index();
  }
}

DEFINE_ASSEMBLER_PLUGIN(
    SumUpElementIndicesPluginNonConst, sumUpElementIndicesNonConst,
    (typename SizeTypeTrait<Assembler>::type *const output)) {
  for (const auto &annotatedElement : this->assembler().meshElements()) {
    *output += annotatedElement.meshView().index();
  }
}

template <class Policy> struct Assembler_Test : Test {
  using Element = NiceMock<test::Element>;
  using Plugins = FeaturePlugins<
      plugins::AssembleEnergyPlugin, plugins::AssembleForceVectorPlugin,
      plugins::AssembleStiffnessMatrixPlugin,
      plugins::UpdateInternalVariablesPlugin,
      plugins::AssembleConsistentMassMatrixPlugin,
      plugins::AssembleLumpedMassMatrixPlugin, plugins::AssembleForceIf,
      SumUpElementIndicesPlugin, SumUpElementIndicesPluginNonConst>;
  using AssemblerWithPlugins = Assembler<Element, Plugins, Policy>;

  using mesh_type = cpppetsc::Mesh<Policy>;

  const double time = .7;

  static constexpr typename mesh_type::size_type numberOfNodes = 2;
  const std::vector<std::vector<int>> connectivity = {
      std::vector<int>{0}, std::vector<int>{0}, std::vector<int>{1}};
  const mesh_type mesh =
      mesh_type::fromConnectivity(Element::SpatialDimension, connectivity,
                                  numberOfNodes, Element::DegreesOfFreedom);
  AssemblerWithPlugins assembler;

  explicit Assembler_Test() {
    for (const auto &element : mesh.localElements()) {
      assembler.emplaceElement(element);
    }
  }

  /**
   * @brief Returns an example global input vector.
   */
  cpppetsc::distributed<typename mesh_type::vector_type> globalInput() const {
    auto result = mesh_type::vector_type::fromGlobalMesh(mesh);
    result.unwrap().replace().elements({0, 1}, {.7, .77});
    return result;
  }

  /**
   * @brief Returns the local part of the example global input.
   */
  cpppetsc::local<typename mesh_type::vector_type> localInput() const {
    auto result = mesh_type::vector_type::fromLocalMesh(mesh);
    mesh.copyToLocalVector(globalInput(), &result);
    return result;
  }

  /**
   * @brief Returns the example global input in canonical order.
   */
  cpppetsc::global<typename mesh_type::vector_type> fullInput() const {
    return mesh_type::vector_type::fromDistributedInCanonicalOrder(
        globalInput(), mesh);
  }
};

template <class Policy>
constexpr typename Assembler_Test<Policy>::mesh_type::size_type
    Assembler_Test<Policy>::numberOfNodes;

using Policies =
    Types<cpppetsc::SequentialComputePolicy, cpppetsc::ParallelComputePolicy>;
TYPED_TEST_CASE(Assembler_Test, Policies);

MATCHER_P(IsWrappedSingleValue, value,
          std::string(negation ? "not " : "") + "wrapped single value " +
              PrintToString(value)) {
  return arg.size() == 1u && arg[0].size() == 1 &&
         cppptest::almost_equal_to_reference(arg[0](0), value);
}

TYPED_TEST(Assembler_Test, type_traits_are_correct) {
  static_assert(
      std::is_same<
          TypeParam,
          typename PolicyTypeTrait<
              typename TestFixture::AssemblerWithPlugins>::type>::value,
      "Policy_type extracts the policy type.");

  static_assert(
      std::is_same<
          typename TestFixture::mesh_type,
          typename MeshTypeTrait<
              typename TestFixture::AssemblerWithPlugins>::type>::value,
      "Mesh_type extracts the mesh type.");

  static_assert(
      std::is_same<
          typename cpppetsc::Vector<TypeParam>,
          typename VectorTypeTrait<
              typename TestFixture::AssemblerWithPlugins>::type>::value,
      "Vector_type extracts the vector type.");

  static_assert(
      std::is_same<
          typename cpppetsc::Matrix<TypeParam>,
          typename MatrixTypeTrait<
              typename TestFixture::AssemblerWithPlugins>::type>::value,
      "Matrix_type extracts the matrix type.");

  static_assert(
      std::is_same<
          typename cpppetsc::Vector<TypeParam>::size_type,
          typename SizeTypeTrait<
              typename TestFixture::AssemblerWithPlugins>::type>::value,
      "Size_type extracts the size type.");

  static_assert(
      std::is_same<
          typename cpppetsc::Vector<TypeParam>::value_type,
          typename ValueTypeTrait<
              typename TestFixture::AssemblerWithPlugins>::type>::value,
      "Value_type extracts the value type.");

  static_assert(
      std::is_same<
          typename TestFixture::Element,
          typename TestFixture::AssemblerWithPlugins::element_type>::value,
      "Element_type extracts the element type.");

  static_assert(
      std::is_same<
          typename TestFixture::Plugins,
          typename PluginTypeTrait<
              typename TestFixture::AssemblerWithPlugins>::type>::value,
      "Plugin_type extracts the Plugin type.");

  static_assert(
      !IsGroupTypeTrait<typename TestFixture::AssemblerWithPlugins>::value,
      "This assembler is not a group.");
}

TYPED_TEST(Assembler_Test, assembling_energy_works) {
  using value_type = typename TestFixture::mesh_type::value_type;

  auto energy = value_type{0.};
  const auto fullInput = this->fullInput();

  for (const auto &annotatedElement : this->assembler.meshElements()) {
    ON_CALL(annotatedElement.instance(),
            computeEnergy(IsWrappedSingleValue(fullInput(
                              annotatedElement.meshView().index() / 2)),
                          this->time))
        .WillByDefault(Return(annotatedElement.meshView().index() + 1.));
  }

  this->assembler.assembleEnergy(this->localInput(), this->time, &energy);
  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPIU_SCALAR, MPIU_SUM,
                            TypeParam::communicator()),
              Eq(0));

  EXPECT_THAT(energy, DoubleEq(1. + 2. + 3.));
}

TYPED_TEST(Assembler_Test, assembling_forces_works) {
  using vector_type = typename TestFixture::mesh_type::vector_type;

  auto forces = vector_type::fromLocalMesh(this->mesh);
  const auto fullInput = this->fullInput();

  for (const auto &annotatedElement : this->assembler.meshElements()) {
    const auto elementForces = typename TestFixture::Element::Forces(
        {typename TestFixture::Element::Vector(
            annotatedElement.meshView().index() + 1.)});
    ON_CALL(annotatedElement.instance(),
            computeForces(IsWrappedSingleValue(fullInput(
                              annotatedElement.meshView().index() / 2)),
                          this->time))
        .WillByDefault(Return(elementForces));
  }

  this->assembler.assembleForceVector(this->localInput(), this->time, &forces);
  auto globalForces = vector_type::fromGlobalMesh(this->mesh);
  this->mesh.addToGlobalVector(forces, &globalForces);

  const auto result =
      vector_type::fromDistributedInCanonicalOrder(globalForces, this->mesh);
  ASSERT_THAT(result.unwrap(), SizeIs(this->numberOfNodes));
  EXPECT_THAT(result(0), DoubleEq(1. + 2.));
  EXPECT_THAT(result(1), DoubleEq(3.));
}

TYPED_TEST(Assembler_Test, assembling_stiffness_matrix_works) {
  using matrix_type = typename TestFixture::mesh_type::matrix_type;

  auto matrix = matrix_type::fromMesh(this->mesh);
  const auto fullInput = this->fullInput();

  for (const auto &annotatedElement : this->assembler.meshElements()) {
    const auto elementMatrix = typename TestFixture::Element::StiffnessMatrix(
        annotatedElement.meshView().index() + 1.);
    ON_CALL(
        annotatedElement.instance(),
        computeStiffnessMatrix(IsWrappedSingleValue(fullInput(
                                   annotatedElement.meshView().index() / 2)),
                               this->time))
        .WillByDefault(Return(elementMatrix));
  }

  this->assembler.assembleStiffnessMatrix(this->localInput(), this->time,
                                          &matrix);
  matrix.finalize();

  ASSERT_THAT(matrix.size(),
              Pair(Eq(this->numberOfNodes), Eq(this->numberOfNodes)));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 0, 1. + 2.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 1, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 0, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 1, 3.));
}

TYPED_TEST(Assembler_Test, assembling_lumped_mass_matrix_works) {
  using matrix_type = typename TestFixture::mesh_type::matrix_type;

  auto matrix = matrix_type::fromMesh(this->mesh);
  const auto fullInput = this->fullInput();

  for (const auto &annotatedElement : this->assembler.meshElements()) {
    const auto elementMatrix = typename TestFixture::Element::StiffnessMatrix(
        annotatedElement.meshView().index() + 1.);
    ON_CALL(annotatedElement.instance(), computeLumpedMassMatrix())
        .WillByDefault(Return(elementMatrix));
  }

  this->assembler.assembleLumpedMassMatrix(&matrix);
  matrix.finalize();

  ASSERT_THAT(matrix.size(),
              Pair(Eq(this->numberOfNodes), Eq(this->numberOfNodes)));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 0, 1. + 2.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 1, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 0, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 1, 3.));
}

TYPED_TEST(Assembler_Test, assembling_consistent_mass_matrix_works) {
  using matrix_type = typename TestFixture::mesh_type::matrix_type;

  auto matrix = matrix_type::fromMesh(this->mesh);
  const auto fullInput = this->fullInput();

  for (const auto &annotatedElement : this->assembler.meshElements()) {
    const auto elementMatrix = typename TestFixture::Element::StiffnessMatrix(
        annotatedElement.meshView().index() + 1.);
    ON_CALL(annotatedElement.instance(), computeConsistentMassMatrix())
        .WillByDefault(Return(elementMatrix));
  }

  this->assembler.assembleConsistentMassMatrix(&matrix);
  matrix.finalize();

  ASSERT_THAT(matrix.size(),
              Pair(Eq(this->numberOfNodes), Eq(this->numberOfNodes)));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 0, 1. + 2.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(0, 1, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 0, 0.));
  EXPECT_THAT(matrix, AlmostEqIfLocal(1, 1, 3.));
}

TYPED_TEST(Assembler_Test, assembling_force_if_condition_holds_works) {
  using value_type = typename TestFixture::mesh_type::value_type;
  using size_type = typename TestFixture::mesh_type::size_type;

  const auto fullInput = this->fullInput();

  for (const auto &annotatedElement : this->assembler.meshElements()) {
    const auto elementForces = typename TestFixture::Element::Forces(
        {typename TestFixture::Element::Vector(
            annotatedElement.meshView().index() + 1.)});
    ON_CALL(annotatedElement.instance(),
            computeForces(IsWrappedSingleValue(fullInput(
                              annotatedElement.meshView().index() / 2)),
                          this->time))
        .WillByDefault(Return(elementForces));
  }

  auto forces = std::vector<value_type>(TestFixture::Element::DegreesOfFreedom);
  this->assembler.assembleForceIf(
      this->localInput(), this->time,
      [&](size_type, const size_type vertexIndex) { return vertexIndex == 0; },
      forces.data());

  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, forces.data(), forces.size(),
                            MPIU_SCALAR, MPIU_SUM, TypeParam::communicator()),
              Eq(0));
  EXPECT_THAT(forces, ElementsAre(DoubleEq(1. + 2.)));
}

TYPED_TEST(Assembler_Test, updating_internal_variables_works) {
  const auto fullInput = this->fullInput();

  for (auto &annotatedElement : this->assembler.meshElements()) {
    EXPECT_CALL(
        annotatedElement.instance(),
        updateInternalVariables(IsWrappedSingleValue(fullInput(
                                    annotatedElement.meshView().index() / 2)),
                                this->time))
        .Times(1);
  }
  this->assembler.updateInternalVariables(this->localInput(), this->time);
}

TYPED_TEST(Assembler_Test, using_const_plugin_works) {
  using size_type = typename TestFixture::mesh_type::size_type;

  auto sum = size_type{0};
  this->assembler.sumUpElementIndices(&sum);
  ASSERT_THAT(MPI_Allreduce(
                  MPI_IN_PLACE, &sum, 1, MPI_INT, MPI_SUM,
                  PolicyTypeTrait<typename TestFixture::AssemblerWithPlugins>::
                      type::communicator()),
              Eq(0));

  EXPECT_THAT(sum, Eq(0 + 1 + 2));
}

TYPED_TEST(Assembler_Test, using_non_const_plugin_works) {
  using size_type = typename TestFixture::mesh_type::size_type;

  auto sum = size_type{0};
  this->assembler.sumUpElementIndicesNonConst(&sum);
  ASSERT_THAT(MPI_Allreduce(
                  MPI_IN_PLACE, &sum, 1, MPI_INT, MPI_SUM,
                  PolicyTypeTrait<typename TestFixture::AssemblerWithPlugins>::
                      type::communicator()),
              Eq(0));

  EXPECT_THAT(sum, Eq(0 + 1 + 2));
}
} // namespace
} // namespace assembly
} // namespace ae108