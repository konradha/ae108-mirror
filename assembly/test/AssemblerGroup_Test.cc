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
#include "ae108/assembly/AssemblerGroup.h"
#include "ae108/assembly/AssemblerTypeTraits.h"
#include "ae108/assembly/FeaturePlugin.h"
#include "ae108/assembly/FeaturePlugins.h"
#include "ae108/assembly/test/Element.h"
#include "ae108/cpppetsc/MeshDataProvider.h"
#include <array>
#include <gmock/gmock.h>

using testing::_;
using testing::DoubleEq;
using testing::Return;
using testing::StrictMock;
using testing::Test;

namespace ae108 {
namespace assembly {
namespace {

DEFINE_CONST_ASSEMBLER_PLUGIN(AddOnePlugin, addOne, (double *const output)) {
  *output += 1.;
}

struct AssemblerGroup_Test : Test {
  using Element = StrictMock<test::Element>;
  using AssemblerWithPlugin = Assembler<Element, FeaturePlugins<AddOnePlugin>>;
  using AssemblerWithoutPlugin = Assembler<Element>;

  using Mesh = AssemblerWithoutPlugin::mesh_type;
  using Connectivity =
      std::array<std::array<AssemblerWithoutPlugin::size_type, 1>, 1>;
  Mesh mesh = Mesh::fromConnectivity<Connectivity>(3, {{{{0}}}}, 1, 0);

  AssemblerGroup<AssemblerWithoutPlugin, AssemblerWithPlugin> group_1;
  AssemblerGroup<AssemblerWithPlugin, AssemblerWithPlugin> group_2;
};

TEST_F(AssemblerGroup_Test, traits_are_inherited_from_members) {
  using Group = AssemblerGroup<AssemblerWithoutPlugin, AssemblerWithPlugin>;

  static_assert(
      std::is_same<PolicyTypeTrait<Group>::type,
                   PolicyTypeTrait<AssemblerWithoutPlugin>::type>::value,
      "Policy type is inherited.");

  static_assert(
      std::is_same<MeshTypeTrait<Group>::type,
                   MeshTypeTrait<AssemblerWithoutPlugin>::type>::value,
      "Mesh type is inherited.");

  static_assert(
      std::is_same<VectorTypeTrait<Group>::type,
                   VectorTypeTrait<AssemblerWithoutPlugin>::type>::value,
      "Vector type is inherited.");

  static_assert(
      std::is_same<MatrixTypeTrait<Group>::type,
                   MatrixTypeTrait<AssemblerWithoutPlugin>::type>::value,
      "Matrix type is inherited.");

  static_assert(
      std::is_same<SizeTypeTrait<Group>::type,
                   SizeTypeTrait<AssemblerWithoutPlugin>::type>::value,
      "Size type is inherited.");

  static_assert(
      std::is_same<ValueTypeTrait<Group>::type,
                   ValueTypeTrait<AssemblerWithoutPlugin>::type>::value,
      "Value type is inherited.");

  static_assert(
      std::is_same<ElementTypeTrait<Group>::type,
                   ElementTypeTrait<AssemblerWithoutPlugin>::type>::value,
      "Element type is inherited.");
}

TEST_F(AssemblerGroup_Test, adds_one_once) {
  double output = 0;
  group_1.addOne(&output);

  EXPECT_THAT(output, DoubleEq(1));
}

TEST_F(AssemblerGroup_Test, adds_one_twice) {
  double output = 0;
  group_2.addOne(&output);

  EXPECT_THAT(output, DoubleEq(2));
}

DEFINE_CONST_ASSEMBLER_PLUGIN(AddEnergyPlugin, addEnergy,
                              (double *const output)) {
  const auto displacements = test::Element::NodalDisplacements();
  for (const auto &annotatedElement : this->assembler().meshElements()) {
    *output += annotatedElement.instance().computeEnergy(displacements, 0.);
  }
}

TEST_F(AssemblerGroup_Test, accessing_elements_in_plugin_works) {
  using EnergyAssembler = Assembler<Element, FeaturePlugins<AddEnergyPlugin>>;
  AssemblerGroup<EnergyAssembler, EnergyAssembler> group;

  const auto values = std::array<double, 2>({.7, .9});

  group.get<0>().emplaceElement(AssemblerWithoutPlugin::ElementView{
      createDataProviderFromMesh(&mesh), 0});
  group.get<1>().emplaceElement(AssemblerWithoutPlugin::ElementView{
      createDataProviderFromMesh(&mesh), 0});

  EXPECT_CALL(group.get<0>().meshElements().begin()->instance(),
              computeEnergy(_, _))
      .WillOnce(Return(values[0]));
  EXPECT_CALL(group.get<1>().meshElements().begin()->instance(),
              computeEnergy(_, _))
      .WillOnce(Return(values[1]));

  double result = 0.;
  group.addEnergy(&result);

  EXPECT_THAT(result, DoubleEq(values[0] + values[1]));
}
} // namespace
} // namespace assembly
} // namespace ae108
