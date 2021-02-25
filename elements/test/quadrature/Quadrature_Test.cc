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

#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/quadrature/integrate.h"
#include <cmath>
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::DoubleNear;
using testing::Not;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace quadrature {
namespace {

template <std::size_t Dimension_, QuadratureType Type_>
struct IntegrableFunction;

template <std::size_t Dimension_>
struct IntegrableFunction<Dimension_, QuadratureType::Cube> {
  double operator()(const std::size_t order,
                    const tensor::Tensor<double, Dimension_> &x) noexcept {
    auto result = double{};
    for (const auto value : x)
      result += std::pow(value + 1., order);
    result /= Dimension_ * std::pow(2., order + Dimension_) / (order + 1.);
    return result;
  }
};

template <std::size_t Dimension_>
struct IntegrableFunction<Dimension_, QuadratureType::Simplex> {
  double operator()(const std::size_t order,
                    const tensor::Tensor<double, Dimension_> &x) noexcept {
    auto result = double{};
    for (const auto value : x)
      result += std::pow(value, order);
    result /= Dimension_ * tgamma(order + 1.) / tgamma(Dimension_ + order + 1.);
    return result;
  }
};

template <class Configuration> struct Quadrature_Test : Test {
  static constexpr std::size_t Dimension = Configuration::Dimension;
  static constexpr std::size_t Order = Configuration::Order;
  static constexpr QuadratureType Type = Configuration::Type;
  using Function = IntegrableFunction<Dimension, Type>;

  /**
   * @brief Returns a test polynomial of the given order evaluated at x. The
   * integral of this polynomial is 1.
   */
  static double
  polynomial_of_order_at(const std::size_t order,
                         const tensor::Tensor<double, Dimension> &x) noexcept {
    return Function{}(order, x);
  }
};

TYPED_TEST_CASE_P(Quadrature_Test);

TYPED_TEST_P(Quadrature_Test,
             integrates_maximum_order_polynomial_correctly) {
  constexpr auto maximum_order = TestFixture::Order;
  const auto f = [&](std::size_t,
                     const tensor::Tensor<double, TestFixture::Dimension> &x) {
    return this->polynomial_of_order_at(maximum_order, x);
  };

  const auto result =
      integrate<Quadrature<TestFixture::Type, TestFixture::Dimension,
                           TestFixture::Order>>(f, 0.);

  EXPECT_THAT(result, DoubleNear(1., 1e-13));
}

TYPED_TEST_P(Quadrature_Test, integrates_maximum_plus_1_order_polynomial_incorrectly) {
  constexpr auto too_high_order = std::size_t{TestFixture::Order + 1};
  const auto f = [&](std::size_t,
                     const tensor::Tensor<double, TestFixture::Dimension> &x) {
    return this->polynomial_of_order_at(too_high_order, x);
  };

  const auto result =
      integrate<Quadrature<TestFixture::Type, TestFixture::Dimension,
                           TestFixture::Order>>(f, 0.);

  EXPECT_THAT(result, Not(DoubleNear(1., 1e-13)));
}

REGISTER_TYPED_TEST_CASE_P(Quadrature_Test,
                           integrates_maximum_order_polynomial_correctly,
                           integrates_maximum_plus_1_order_polynomial_incorrectly);

template <std::size_t Dimension_, std::size_t Order_, QuadratureType Type_>
struct Configuration {
  static constexpr std::size_t Dimension = Dimension_;
  static constexpr std::size_t Order = Order_;
  static constexpr QuadratureType Type = Type_;
};

using Configurations = Types<Configuration<1, 1, QuadratureType::Cube>,
                             Configuration<2, 1, QuadratureType::Cube>,
                             Configuration<3, 1, QuadratureType::Cube>,
                             Configuration<1, 3, QuadratureType::Cube>,
                             Configuration<2, 3, QuadratureType::Cube>,
                             Configuration<3, 3, QuadratureType::Cube>,
                             Configuration<1, 5, QuadratureType::Cube>,
                             Configuration<2, 5, QuadratureType::Cube>,
                             Configuration<3, 5, QuadratureType::Cube>,
                             Configuration<1, 7, QuadratureType::Cube>,
                             Configuration<2, 7, QuadratureType::Cube>,
                             Configuration<3, 7, QuadratureType::Cube>,
                             Configuration<1, 9, QuadratureType::Cube>,
                             Configuration<2, 9, QuadratureType::Cube>,
                             Configuration<3, 9, QuadratureType::Cube>,
                             Configuration<2, 1, QuadratureType::Simplex>,
                             Configuration<2, 2, QuadratureType::Simplex>,
                             Configuration<2, 3, QuadratureType::Simplex>>;
INSTANTIATE_TYPED_TEST_CASE_P(Quadrature_Test, Quadrature_Test, Configurations);

struct Quadrature_1D_Test : Test {
  static constexpr std::size_t Dimension = 1;
};

TEST_F(Quadrature_1D_Test, passes_on_args) {
  const auto f = [](std::size_t, const tensor::Tensor<double, Dimension> &,
                    const int &arg) { return arg; };

  const Quadrature<QuadratureType::Cube, Dimension, 1>::Collection<int> args = {
      {7}};

  const auto result =
      integrate<Quadrature<QuadratureType::Cube, Dimension, 1>>(f, 0., args);

  EXPECT_THAT(result, DoubleEq(7. * 2.));
}

TEST_F(Quadrature_1D_Test, passes_on_ids) {
  const auto f = [](std::size_t id, const tensor::Tensor<double, Dimension> &) {
    return id + std::size_t{1};
  };

  const auto result =
      integrate<Quadrature<QuadratureType::Cube, Dimension, 3>>(f, 0.);

  EXPECT_THAT(result, DoubleEq((1. + 0.) + (1. + 1.)));
}

TEST_F(Quadrature_1D_Test, adds_to_initial_value) {
  const auto f = [](std::size_t, const tensor::Tensor<double, Dimension> &) {
    return 1.;
  };

  const auto result =
      integrate<Quadrature<QuadratureType::Cube, Dimension, 3>>(f, 7.);

  EXPECT_THAT(result, DoubleEq(2. * 1. + 7.));
}

} // namespace
} // namespace quadrature
} // namespace elements
} // namespace ae108