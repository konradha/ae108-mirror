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
#include "ae108/elements/materialmodels/Minimal.h"
#include "ae108/elements/materialmodels/compute_energy.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace materialmodels {
namespace {

template <class MaterialModel_> struct Configuration {
  using MaterialModel = MaterialModel_;
  static MaterialModel create_model() noexcept { return MaterialModel(); }

  static typename MaterialModel_::Time create_time() noexcept {
    return typename MaterialModel::Time{0.};
  }
};

using Configurations =
    Types<Configuration<Minimal<1>>, Configuration<Minimal<2>>,
          Configuration<Minimal<3>>>;
INSTANTIATE_TYPED_TEST_CASE_P(Minimal_Test, MaterialModel_Test, Configurations);

struct Minimal_1D_Test : Test {
  using MaterialModel = Minimal<1>;
  MaterialModel model = MaterialModel();
};

TEST_F(Minimal_1D_Test, computes_zero_energy) {
  const auto time = MaterialModel::Time{0};
  const auto gradient = MaterialModel::DisplacementGradient();

  const auto result = compute_energy(model, 0, gradient, time);

  EXPECT_THAT(result, DoubleEq(0.));
}

} // namespace
} // namespace materialmodels
} // namespace elements
} // namespace ae108