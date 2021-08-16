// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "Element_Test.h"
#include "ae108/elements/CoreElement.h"
#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/integrator/IsoparametricIntegrator.h"
#include "ae108/elements/materialmodels/AutomaticStressTrait.h"
#include "ae108/elements/materialmodels/AutomaticTangentMatrixTrait.h"
#include "ae108/elements/materialmodels/Hookean.h"
#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/shape/Hexa8.h"
#include "ae108/elements/shape/Quad4.h"
#include "ae108/elements/shape/Seg2.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace {

template <class Element_> struct Configuration_1D {
  using Element = Element_;
  static Element create_element() noexcept {
    using Integrator = typename Element::Integrator;
    using MaterialModel = typename Element::MaterialModel;
    using Embedding = typename Integrator::Embedding;

    return Element(MaterialModel(1., 0.), Integrator(Embedding({{
                                              {{0.}},
                                              {{1.}},
                                          }})));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Element_Seg2 = CoreElement<
    materialmodels::Hookean<1>,
    integrator::IsoparametricIntegrator<
        shape::Seg2,
        quadrature::Quadrature<quadrature::QuadratureType::Cube, 1, 1>>>;

using Configurations_1D = Types<Configuration_1D<Element_Seg2>>;
INSTANTIATE_TYPED_TEST_CASE_P(CoreElement_Seg2_Test, Element_Test,
                              Configurations_1D);

struct CoreElement_Seg2_Test : Test {
  using Element = Element_Seg2;
  const Element element = Configuration_1D<Element>::create_element();
};

TEST_F(CoreElement_Seg2_Test, computes_energy_with_no_displacements) {
  const auto time = Element::Time{0.};
  const auto displacements = Element::NodalDisplacements();

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(0.));
}

TEST_F(CoreElement_Seg2_Test, computes_energy_with_displacements_1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0.}},
      {{1.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(.5));
}

TEST_F(CoreElement_Seg2_Test, computes_energy_with_displacements_2) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0.}},
      {{2.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(2.));
}

template <class Element_> struct Configuration_2D {
  using Element = Element_;
  static Element create_element() noexcept {
    using Integrator = typename Element::Integrator;
    using MaterialModel = typename Element::MaterialModel;
    using Embedding = typename Integrator::Embedding;

    return Element(MaterialModel(1., 0.), Integrator(Embedding({{
                                              {{0., 0.}},
                                              {{1., 0.}},
                                              {{1., 1.}},
                                              {{0., 1.}},
                                          }})));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Element_Quad4 = CoreElement<
    materialmodels::Hookean<2>,
    integrator::IsoparametricIntegrator<
        shape::Quad4,
        quadrature::Quadrature<quadrature::QuadratureType::Cube, 2, 1>>>;

using Configurations_2D = Types<Configuration_2D<Element_Quad4>>;
INSTANTIATE_TYPED_TEST_CASE_P(CoreElement_Quad4_Test, Element_Test,
                              Configurations_2D);

struct CoreElement_Quad4_Test : Test {
  using Element = Element_Quad4;
  const Element element = Configuration_2D<Element>::create_element();
};

TEST_F(CoreElement_Quad4_Test, computes_energy_with_no_displacements) {
  const auto time = Element::Time{0.};
  const auto displacements = Element::NodalDisplacements();

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(0.));
}

TEST_F(CoreElement_Quad4_Test, computes_correct_energy_with_displacements_1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{1., 0.}},
      {{1., 0.}},
      {{0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(.5));
}

TEST_F(CoreElement_Quad4_Test, computes_correct_energy_with_displacements_2) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{-2., 0.}},
      {{+0., 0.}},
      {{+0., 0.}},
      {{-2., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(2.));
}

template <class Element_> struct Configuration_3D {
  using Element = Element_;
  static Element create_element() noexcept {
    using Integrator = typename Element::Integrator;
    using MaterialModel = typename Element::MaterialModel;
    using Embedding = typename Integrator::Embedding;

    return Element(MaterialModel(1., 0.), Integrator(Embedding({{
                                              {{0., 0., 0.}},
                                              {{1., 0., 0.}},
                                              {{1., 1., 0.}},
                                              {{0., 1., 0.}},
                                              {{0., 0., 1.}},
                                              {{1., 0., 1.}},
                                              {{1., 1., 1.}},
                                              {{0., 1., 1.}},
                                          }})));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Element_Hexa8 = CoreElement<
    materialmodels::Hookean<3>,
    integrator::IsoparametricIntegrator<
        shape::Hexa8,
        quadrature::Quadrature<quadrature::QuadratureType::Cube, 3, 1>>>;

using Configurations_3D = Types<Configuration_3D<Element_Hexa8>>;
INSTANTIATE_TYPED_TEST_CASE_P(CoreElement_Hexa8_Test, Element_Test,
                              Configurations_3D);

struct CoreElement_Hexa8_Test : Test {
  using Element = Element_Hexa8;
  const Element element = Configuration_3D<Element>::create_element();
};

TEST_F(CoreElement_Hexa8_Test, computes_correct_energy_with_displacements_1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0., 0.}},
      {{1., 0., 0.}},
      {{1., 0., 0.}},
      {{0., 0., 0.}},
      {{0., 0., 0.}},
      {{1., 0., 0.}},
      {{1., 0., 0.}},
      {{0., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(.5));
}

TEST_F(CoreElement_Hexa8_Test, computes_correct_energy_with_displacements_2) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{2., 0., 0.}},
      {{0., 0., 0.}},
      {{0., 0., 0.}},
      {{2., 0., 0.}},
      {{2., 0., 0.}},
      {{0., 0., 0.}},
      {{0., 0., 0.}},
      {{2., 0., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(2.));
}

/**
 * @brief A material model with a number of degrees of freedom different from
 * the dimension.
 */
struct Model_Dof final
    : materialmodels::MaterialModelBase<std::size_t, double, 1 /* dimension */,
                                        2 /* degrees of freedom */> {
  template <class... Args> explicit Model_Dof(Args &&...) {}
};
} // namespace

namespace materialmodels {

/**
 * @brief Returns (x + y / 2)^2.
 */
template <> struct ComputeEnergyTrait<Model_Dof> {
  template <class MaterialModel>
  typename MaterialModel::Energy
  operator()(const MaterialModel &, const typename MaterialModel::size_type,
             const typename MaterialModel::DisplacementGradient &x,
             const typename MaterialModel::Time) noexcept {
    return std::pow(x.at(0).at(0) + x.at(1).at(0) / 2., 2.);
  }
};

template <>
struct ComputeStressTrait<Model_Dof> : AutomaticStressTrait<Model_Dof> {};

template <>
struct ComputeTangentMatrixTrait<Model_Dof>
    : AutomaticTangentMatrixTrait<Model_Dof> {};

} // namespace materialmodels

namespace {
using Element_Dof = CoreElement<
    Model_Dof,
    integrator::IsoparametricIntegrator<
        shape::Seg2,
        quadrature::Quadrature<quadrature::QuadratureType::Cube, 1, 1>>>;

using Configurations_Dof = Types<Configuration_1D<Element_Dof>>;
INSTANTIATE_TYPED_TEST_CASE_P(CoreElement_Dof_Test, Element_Test,
                              Configurations_Dof);

struct CoreElement_Dof_Test : Test {
  using Element = Element_Dof;
  const Element element = Configuration_1D<Element>::create_element();
};

TEST_F(CoreElement_Dof_Test, computes_energy_with_x_displacements_1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{1., 0.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(1.));
}

TEST_F(CoreElement_Dof_Test, computes_energy_with_y_displacements_1) {
  const auto time = Element::Time{0.};

  const Element::NodalDisplacements displacements = {{
      {{0., 0.}},
      {{0., 1.}},
  }};

  EXPECT_THAT(element.computeEnergy(displacements, time), DoubleEq(.25));
}

} // namespace
} // namespace elements
} // namespace ae108