// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/integrator/IsoparametricIntegrator.h"
#include "ae108/elements/integrator/compute_volume.h"
#include "ae108/elements/integrator/integrate.h"
#include "ae108/elements/integrator/integrate_shape.h"
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