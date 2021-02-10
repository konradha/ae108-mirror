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

#include "MaterialModel_Test.h"
#include "ae108/elements/materialmodels/Hookean.h"
#include "ae108/elements/materialmodels/compute_energy.h"
#include "ae108/elements/materialmodels/compute_strain.h"
#include <gmock/gmock.h>
#include <type_traits>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace materialmodels {
namespace {

template <class MaterialModel_> struct Configuration {
  using MaterialModel = MaterialModel_;
  static MaterialModel create_model() noexcept {
    constexpr auto youngs_modulus = 3.;
    constexpr auto poission_ratio = .2;
    return MaterialModel(youngs_modulus, poission_ratio);
  }

  static typename MaterialModel_::Time create_time() noexcept {
    return typename MaterialModel::Time{0.};
  }
};

using Configurations =
    Types<Configuration<Hookean<1>>, Configuration<Hookean<2>>,
          Configuration<Hookean<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(Hookean_Test, MaterialModel_Test, Configurations);

struct Hookean_1D_Test : Test {
  const double youngs_modulus = 3.;
  const double poisson_ratio = .2;
  using MaterialModel = Hookean<1>;
  MaterialModel model = MaterialModel(youngs_modulus, poisson_ratio);
};

TEST_F(Hookean_1D_Test, correct_lambda) {
  EXPECT_THAT(model.lambda(), DoubleEq(5. / 6.));
}

TEST_F(Hookean_1D_Test, correct_mu) {
  EXPECT_THAT(model.mu(), DoubleEq(5. / 4.));
}

TEST_F(Hookean_1D_Test, zero_energy_for_no_gradient) {
  const auto time = MaterialModel::Time{0};
  const auto gradient = MaterialModel::DisplacementGradient();

  const auto result =
      compute_energy(model, MaterialModel::unknown_id(), gradient, time);

  EXPECT_THAT(result, DoubleEq(0.));
}

TEST_F(Hookean_1D_Test, zero_strain_for_no_gradient) {
  const auto time = MaterialModel::Time{0};
  const auto gradient = MaterialModel::DisplacementGradient();

  const auto result =
      compute_strain(model, MaterialModel::unknown_id(), gradient, time);

  EXPECT_THAT(result, ElementsAre(ElementsAre(DoubleEq(0.))));
}

TEST_F(Hookean_1D_Test, correct_energy_for_gradient_1) {
  const auto time = MaterialModel::Time{0};
  const MaterialModel::DisplacementGradient gradient = {{
      {{1.}},
  }};

  const auto result =
      compute_energy(model, MaterialModel::unknown_id(), gradient, time);

  EXPECT_THAT(result, DoubleEq(model.mu() + .5 * model.lambda()));
}

TEST_F(Hookean_1D_Test, correct_strain_for_gradient_1) {
  const auto time = MaterialModel::Time{0};
  const MaterialModel::DisplacementGradient gradient = {{
      {{1.}},
  }};

  const auto result =
      compute_strain(model, MaterialModel::unknown_id(), gradient, time);

  EXPECT_THAT(result, ElementsAre(ElementsAre(DoubleEq(1.))));
}

TEST_F(Hookean_1D_Test, correct_energy_for_gradient_2) {
  const auto time = MaterialModel::Time{0};
  const MaterialModel::DisplacementGradient gradient = {{
      {{2.}},
  }};

  const auto result =
      compute_energy(model, MaterialModel::unknown_id(), gradient, time);

  EXPECT_THAT(result, DoubleEq(2. * (2. * model.mu() + model.lambda())));
}

TEST_F(Hookean_1D_Test, correct_strain_for_gradient_2) {
  const auto time = MaterialModel::Time{0};
  const MaterialModel::DisplacementGradient gradient = {{
      {{2.}},
  }};

  const auto result = compute_strain(model, 0, gradient, time);

  EXPECT_THAT(result, ElementsAre(ElementsAre(DoubleEq(2.))));
}

struct Hookean_3D_Test : Test {
  const double youngs_modulus = 3.;
  const double poisson_ratio = .2;
  using MaterialModel = Hookean<3>;
  MaterialModel model = MaterialModel(youngs_modulus, poisson_ratio);
};

TEST_F(Hookean_3D_Test, symmetrized_gradient_is_strain) {
  const auto time = MaterialModel::Time{0};
  const auto gradient = MaterialModel::DisplacementGradient({{
      {{1., 0., 3.}},
      {{4., 3., 8.}},
      {{3., 0., 5.}},
  }});
  const auto result =
      compute_strain(model, MaterialModel::unknown_id(), gradient, time);

  EXPECT_THAT(
      result,
      ElementsAre(ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.)),
                  ElementsAre(DoubleEq(2.), DoubleEq(3.), DoubleEq(4.)),
                  ElementsAre(DoubleEq(3.), DoubleEq(4.), DoubleEq(5.))));
}

TEST_F(Hookean_3D_Test, correct_energy_for_identity_gradient) {
  const auto time = MaterialModel::Time{0};
  const auto gradient = MaterialModel::DisplacementGradient({{
      {{1., 0., 0.}},
      {{0., 1., 0.}},
      {{0., 0., 1.}},
  }});
  const auto result = compute_energy(model, 0, gradient, time);

  EXPECT_THAT(result, DoubleEq(model.mu() * 3. + .5 * model.lambda() * 9.));
}

TEST_F(Hookean_3D_Test, correct_energy_for_upper_diagonal_gradient) {
  const auto time = MaterialModel::Time{0};
  const auto gradient = MaterialModel::DisplacementGradient({{
      {{1., 1., 1.}},
      {{0., 1., 1.}},
      {{0., 0., 1.}},
  }});
  const auto result = compute_energy(model, 0, gradient, time);

  EXPECT_THAT(result, DoubleEq(model.mu() * (3. + .25 * 6.) +
                               .5 * model.lambda() * 9.));
}

} // namespace
} // namespace materialmodels
} // namespace elements
} // namespace ae108