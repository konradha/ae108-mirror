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