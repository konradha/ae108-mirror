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

#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/integrator/IsoparametricIntegrator.h"
#include "ae108/elements/integrator/integrate.h"
#include "ae108/elements/integrator/integrate_shape.h"
#include "ae108/elements/integrator/volume.h"
#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/shape/Seg2.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace elements {
namespace integrator {
namespace {

struct IsoperimetricIntegrator_Test : Test {
  using Shape = shape::Seg2;
  using Embedding = embedding::IsoparametricEmbedding<Shape>;
  using Quadrature =
      quadrature::Quadrature<quadrature::QuadratureType::Cube, 1, 1>;
  using Integrator = IsoparametricIntegrator<Shape, Quadrature>;

  const Embedding embedding = Embedding{{{
      {{1.}},
      {{2.}},
  }}};
  const Integrator integrator = Integrator(embedding);
};

TEST_F(IsoperimetricIntegrator_Test, pre_transform_is_correct) {
  const auto &pre = integrator.pre();

  EXPECT_THAT(pre, ElementsAre(ElementsAre(ElementsAre(DoubleEq(-1.)),
                                           ElementsAre(DoubleEq(1.)))));
}

TEST_F(IsoperimetricIntegrator_Test, post_transform_is_correct) {
  const auto &post = integrator.post();

  EXPECT_THAT(post, ElementsAre(DoubleEq(.5)));
}

TEST_F(IsoperimetricIntegrator_Test, computing_volume_via_integration_works) {
  const auto result = integrate(
      integrator,
      [](const Integrator::size_type &, const Integrator::Point<1> &,
         const Integrator::PreTransform &) { return 1.; },
      Integrator::DiscretizedFunction<1>(), 0.);

  EXPECT_THAT(result, DoubleEq(1.));
}

TEST_F(IsoperimetricIntegrator_Test, computing_volume_works) {
  const auto result = volume(integrator);
  EXPECT_THAT(result, DoubleEq(1.));
}

TEST_F(IsoperimetricIntegrator_Test, integrating_shape_works) {
  const auto factor = std::sqrt(2.);
  const auto result = integrate_shape(
      integrator,
      [factor](auto &&, const auto &point) {
        return point[0] + factor * point[1];
      },
      0.);

  EXPECT_THAT(result, DoubleEq(1. / 2. + factor * 1. / 2.));
}

} // namespace
} // namespace integrator
} // namespace elements
} // namespace ae108