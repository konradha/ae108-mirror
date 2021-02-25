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

#pragma once

#include "ae108/elements/quadrature/IntegrateTrait.h"
#include "ae108/elements/quadrature/QuadratureBase.h"
#include "ae108/elements/tensor/Tensor.h"
#include <cmath>
#include <cstddef>

namespace ae108 {
namespace elements {
namespace quadrature {

enum class QuadratureType { Cube, Simplex };

/**
 * @brief Contains a static member "data" that contains the points and weights.
 * @tparam Order_ The maximum order of a polynomial that the quadrature rule
 * integrates exactly.
 * @remark Derives from QuadratureBase.
 */
template <QuadratureType Type_, std::size_t Dimension_, std::size_t Order_>
struct Quadrature;

#define AE108_ELEMENTS_QUADRATURE_DEFINE(type, dimension, order, size, ...)    \
  template <>                                                                  \
  struct Quadrature<type, dimension, order>                                    \
      : QuadratureBase<std::size_t, tensor::Tensor<double, dimension>, double, \
                       size> {                                                 \
    static constexpr Data data = __VA_ARGS__;                                  \
  }

AE108_ELEMENTS_QUADRATURE_DEFINE(QuadratureType::Cube, 1, 1, 1,
                                 {{{
                                      {{+0.0000000000000000}},
                                  }},
                                  {{+2.0000000000000000}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 2, 1, 1,
    {{{
         {{+0.0000000000000000, +0.0000000000000000}},
     }},
     {{+4.0000000000000000}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 3, 1, 1,
    {{{
         {{+0.0000000000000000, +0.0000000000000000, +0.0000000000000000}},
     }},
     {{+8.0000000000000000}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(QuadratureType::Cube, 1, 3, 2,
                                 {{{
                                      {{-0.5773502691896257}},
                                      {{+0.5773502691896257}},
                                  }},
                                  {{+1.0000000000000000,
                                    +1.0000000000000000}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 2, 3, 4,
    {{{
         {{-0.5773502691896257, -0.5773502691896257}},
         {{-0.5773502691896257, +0.5773502691896257}},
         {{+0.5773502691896257, -0.5773502691896257}},
         {{+0.5773502691896257, +0.5773502691896257}},
     }},
     {{+1.0000000000000000, +1.0000000000000000, +1.0000000000000000,
       +1.0000000000000000}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 3, 3, 8,
    {{{
         {{-0.5773502691896257, -0.5773502691896257, -0.5773502691896257}},
         {{-0.5773502691896257, -0.5773502691896257, +0.5773502691896257}},
         {{-0.5773502691896257, +0.5773502691896257, -0.5773502691896257}},
         {{-0.5773502691896257, +0.5773502691896257, +0.5773502691896257}},
         {{+0.5773502691896257, -0.5773502691896257, -0.5773502691896257}},
         {{+0.5773502691896257, -0.5773502691896257, +0.5773502691896257}},
         {{+0.5773502691896257, +0.5773502691896257, -0.5773502691896257}},
         {{+0.5773502691896257, +0.5773502691896257, +0.5773502691896257}},
     }},
     {{+1.0000000000000000, +1.0000000000000000, +1.0000000000000000,
       +1.0000000000000000, +1.0000000000000000, +1.0000000000000000,
       +1.0000000000000000, +1.0000000000000000}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(QuadratureType::Cube, 1, 5, 3,
                                 {{{
                                      {{-0.7745966692414834}},
                                      {{+0.0000000000000000}},
                                      {{+0.7745966692414834}},
                                  }},
                                  {{+0.5555555555555556, +0.8888888888888888,
                                    +0.5555555555555556}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 2, 5, 9,
    {{{
         {{-0.7745966692414834, -0.7745966692414834}},
         {{-0.7745966692414834, +0.0000000000000000}},
         {{-0.7745966692414834, +0.7745966692414834}},
         {{+0.0000000000000000, -0.7745966692414834}},
         {{+0.0000000000000000, +0.0000000000000000}},
         {{+0.0000000000000000, +0.7745966692414834}},
         {{+0.7745966692414834, -0.7745966692414834}},
         {{+0.7745966692414834, +0.0000000000000000}},
         {{+0.7745966692414834, +0.7745966692414834}},
     }},
     {{+0.3086419753086420, +0.4938271604938271, +0.3086419753086420,
       +0.4938271604938271, +0.7901234567901234, +0.4938271604938271,
       +0.3086419753086420, +0.4938271604938271, +0.3086419753086420}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 3, 5, 27,
    {{{
         {{-0.7745966692414834, -0.7745966692414834, -0.7745966692414834}},
         {{-0.7745966692414834, -0.7745966692414834, +0.0000000000000000}},
         {{-0.7745966692414834, -0.7745966692414834, +0.7745966692414834}},
         {{-0.7745966692414834, +0.0000000000000000, -0.7745966692414834}},
         {{-0.7745966692414834, +0.0000000000000000, +0.0000000000000000}},
         {{-0.7745966692414834, +0.0000000000000000, +0.7745966692414834}},
         {{-0.7745966692414834, +0.7745966692414834, -0.7745966692414834}},
         {{-0.7745966692414834, +0.7745966692414834, +0.0000000000000000}},
         {{-0.7745966692414834, +0.7745966692414834, +0.7745966692414834}},
         {{+0.0000000000000000, -0.7745966692414834, -0.7745966692414834}},
         {{+0.0000000000000000, -0.7745966692414834, +0.0000000000000000}},
         {{+0.0000000000000000, -0.7745966692414834, +0.7745966692414834}},
         {{+0.0000000000000000, +0.0000000000000000, -0.7745966692414834}},
         {{+0.0000000000000000, +0.0000000000000000, +0.0000000000000000}},
         {{+0.0000000000000000, +0.0000000000000000, +0.7745966692414834}},
         {{+0.0000000000000000, +0.7745966692414834, -0.7745966692414834}},
         {{+0.0000000000000000, +0.7745966692414834, +0.0000000000000000}},
         {{+0.0000000000000000, +0.7745966692414834, +0.7745966692414834}},
         {{+0.7745966692414834, -0.7745966692414834, -0.7745966692414834}},
         {{+0.7745966692414834, -0.7745966692414834, +0.0000000000000000}},
         {{+0.7745966692414834, -0.7745966692414834, +0.7745966692414834}},
         {{+0.7745966692414834, +0.0000000000000000, -0.7745966692414834}},
         {{+0.7745966692414834, +0.0000000000000000, +0.0000000000000000}},
         {{+0.7745966692414834, +0.0000000000000000, +0.7745966692414834}},
         {{+0.7745966692414834, +0.7745966692414834, -0.7745966692414834}},
         {{+0.7745966692414834, +0.7745966692414834, +0.0000000000000000}},
         {{+0.7745966692414834, +0.7745966692414834, +0.7745966692414834}},
     }},
     {{+0.1714677640603567, +0.2743484224965707, +0.1714677640603567,
       +0.2743484224965707, +0.4389574759945130, +0.2743484224965707,
       +0.1714677640603567, +0.2743484224965707, +0.1714677640603567,
       +0.2743484224965707, +0.4389574759945130, +0.2743484224965707,
       +0.4389574759945130, +0.7023319615912208, +0.4389574759945130,
       +0.2743484224965707, +0.4389574759945130, +0.2743484224965707,
       +0.1714677640603567, +0.2743484224965707, +0.1714677640603567,
       +0.2743484224965707, +0.4389574759945130, +0.2743484224965707,
       +0.1714677640603567, +0.2743484224965707, +0.1714677640603567}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(QuadratureType::Cube, 1, 7, 4,
                                 {{{
                                      {{-0.8611363115940526}},
                                      {{-0.3399810435848563}},
                                      {{+0.3399810435848563}},
                                      {{+0.8611363115940526}},
                                  }},
                                  {{+0.3478548451374538, +0.6521451548625461,
                                    +0.6521451548625461,
                                    +0.3478548451374538}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 2, 7, 16,
    {{{
         {{-0.8611363115940526, -0.8611363115940526}},
         {{-0.8611363115940526, -0.3399810435848563}},
         {{-0.8611363115940526, +0.3399810435848563}},
         {{-0.8611363115940526, +0.8611363115940526}},
         {{-0.3399810435848563, -0.8611363115940526}},
         {{-0.3399810435848563, -0.3399810435848563}},
         {{-0.3399810435848563, +0.3399810435848563}},
         {{-0.3399810435848563, +0.8611363115940526}},
         {{+0.3399810435848563, -0.8611363115940526}},
         {{+0.3399810435848563, -0.3399810435848563}},
         {{+0.3399810435848563, +0.3399810435848563}},
         {{+0.3399810435848563, +0.8611363115940526}},
         {{+0.8611363115940526, -0.8611363115940526}},
         {{+0.8611363115940526, -0.3399810435848563}},
         {{+0.8611363115940526, +0.3399810435848563}},
         {{+0.8611363115940526, +0.8611363115940526}},
     }},
     {{+0.1210029932856020, +0.2268518518518519, +0.2268518518518519,
       +0.1210029932856020, +0.2268518518518519, +0.4252933030106943,
       +0.4252933030106943, +0.2268518518518519, +0.2268518518518519,
       +0.4252933030106943, +0.4252933030106943, +0.2268518518518519,
       +0.1210029932856020, +0.2268518518518519, +0.2268518518518519,
       +0.1210029932856020}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 3, 7, 64,
    {{{
         {{-0.8611363115940526, -0.8611363115940526, -0.8611363115940526}},
         {{-0.8611363115940526, -0.8611363115940526, -0.3399810435848563}},
         {{-0.8611363115940526, -0.8611363115940526, +0.3399810435848563}},
         {{-0.8611363115940526, -0.8611363115940526, +0.8611363115940526}},
         {{-0.8611363115940526, -0.3399810435848563, -0.8611363115940526}},
         {{-0.8611363115940526, -0.3399810435848563, -0.3399810435848563}},
         {{-0.8611363115940526, -0.3399810435848563, +0.3399810435848563}},
         {{-0.8611363115940526, -0.3399810435848563, +0.8611363115940526}},
         {{-0.8611363115940526, +0.3399810435848563, -0.8611363115940526}},
         {{-0.8611363115940526, +0.3399810435848563, -0.3399810435848563}},
         {{-0.8611363115940526, +0.3399810435848563, +0.3399810435848563}},
         {{-0.8611363115940526, +0.3399810435848563, +0.8611363115940526}},
         {{-0.8611363115940526, +0.8611363115940526, -0.8611363115940526}},
         {{-0.8611363115940526, +0.8611363115940526, -0.3399810435848563}},
         {{-0.8611363115940526, +0.8611363115940526, +0.3399810435848563}},
         {{-0.8611363115940526, +0.8611363115940526, +0.8611363115940526}},
         {{-0.3399810435848563, -0.8611363115940526, -0.8611363115940526}},
         {{-0.3399810435848563, -0.8611363115940526, -0.3399810435848563}},
         {{-0.3399810435848563, -0.8611363115940526, +0.3399810435848563}},
         {{-0.3399810435848563, -0.8611363115940526, +0.8611363115940526}},
         {{-0.3399810435848563, -0.3399810435848563, -0.8611363115940526}},
         {{-0.3399810435848563, -0.3399810435848563, -0.3399810435848563}},
         {{-0.3399810435848563, -0.3399810435848563, +0.3399810435848563}},
         {{-0.3399810435848563, -0.3399810435848563, +0.8611363115940526}},
         {{-0.3399810435848563, +0.3399810435848563, -0.8611363115940526}},
         {{-0.3399810435848563, +0.3399810435848563, -0.3399810435848563}},
         {{-0.3399810435848563, +0.3399810435848563, +0.3399810435848563}},
         {{-0.3399810435848563, +0.3399810435848563, +0.8611363115940526}},
         {{-0.3399810435848563, +0.8611363115940526, -0.8611363115940526}},
         {{-0.3399810435848563, +0.8611363115940526, -0.3399810435848563}},
         {{-0.3399810435848563, +0.8611363115940526, +0.3399810435848563}},
         {{-0.3399810435848563, +0.8611363115940526, +0.8611363115940526}},
         {{+0.3399810435848563, -0.8611363115940526, -0.8611363115940526}},
         {{+0.3399810435848563, -0.8611363115940526, -0.3399810435848563}},
         {{+0.3399810435848563, -0.8611363115940526, +0.3399810435848563}},
         {{+0.3399810435848563, -0.8611363115940526, +0.8611363115940526}},
         {{+0.3399810435848563, -0.3399810435848563, -0.8611363115940526}},
         {{+0.3399810435848563, -0.3399810435848563, -0.3399810435848563}},
         {{+0.3399810435848563, -0.3399810435848563, +0.3399810435848563}},
         {{+0.3399810435848563, -0.3399810435848563, +0.8611363115940526}},
         {{+0.3399810435848563, +0.3399810435848563, -0.8611363115940526}},
         {{+0.3399810435848563, +0.3399810435848563, -0.3399810435848563}},
         {{+0.3399810435848563, +0.3399810435848563, +0.3399810435848563}},
         {{+0.3399810435848563, +0.3399810435848563, +0.8611363115940526}},
         {{+0.3399810435848563, +0.8611363115940526, -0.8611363115940526}},
         {{+0.3399810435848563, +0.8611363115940526, -0.3399810435848563}},
         {{+0.3399810435848563, +0.8611363115940526, +0.3399810435848563}},
         {{+0.3399810435848563, +0.8611363115940526, +0.8611363115940526}},
         {{+0.8611363115940526, -0.8611363115940526, -0.8611363115940526}},
         {{+0.8611363115940526, -0.8611363115940526, -0.3399810435848563}},
         {{+0.8611363115940526, -0.8611363115940526, +0.3399810435848563}},
         {{+0.8611363115940526, -0.8611363115940526, +0.8611363115940526}},
         {{+0.8611363115940526, -0.3399810435848563, -0.8611363115940526}},
         {{+0.8611363115940526, -0.3399810435848563, -0.3399810435848563}},
         {{+0.8611363115940526, -0.3399810435848563, +0.3399810435848563}},
         {{+0.8611363115940526, -0.3399810435848563, +0.8611363115940526}},
         {{+0.8611363115940526, +0.3399810435848563, -0.8611363115940526}},
         {{+0.8611363115940526, +0.3399810435848563, -0.3399810435848563}},
         {{+0.8611363115940526, +0.3399810435848563, +0.3399810435848563}},
         {{+0.8611363115940526, +0.3399810435848563, +0.8611363115940526}},
         {{+0.8611363115940526, +0.8611363115940526, -0.8611363115940526}},
         {{+0.8611363115940526, +0.8611363115940526, -0.3399810435848563}},
         {{+0.8611363115940526, +0.8611363115940526, +0.3399810435848563}},
         {{+0.8611363115940526, +0.8611363115940526, +0.8611363115940526}},
     }},
     {{+0.0420914774905315, +0.0789115157950706, +0.0789115157950706,
       +0.0420914774905315, +0.0789115157950706, +0.1479403360567813,
       +0.1479403360567813, +0.0789115157950706, +0.0789115157950706,
       +0.1479403360567813, +0.1479403360567813, +0.0789115157950706,
       +0.0420914774905315, +0.0789115157950706, +0.0789115157950706,
       +0.0420914774905315, +0.0789115157950706, +0.1479403360567813,
       +0.1479403360567813, +0.0789115157950706, +0.1479403360567813,
       +0.2773529669539130, +0.2773529669539130, +0.1479403360567813,
       +0.1479403360567813, +0.2773529669539130, +0.2773529669539130,
       +0.1479403360567813, +0.0789115157950706, +0.1479403360567813,
       +0.1479403360567813, +0.0789115157950706, +0.0789115157950706,
       +0.1479403360567813, +0.1479403360567813, +0.0789115157950706,
       +0.1479403360567813, +0.2773529669539130, +0.2773529669539130,
       +0.1479403360567813, +0.1479403360567813, +0.2773529669539130,
       +0.2773529669539130, +0.1479403360567813, +0.0789115157950706,
       +0.1479403360567813, +0.1479403360567813, +0.0789115157950706,
       +0.0420914774905315, +0.0789115157950706, +0.0789115157950706,
       +0.0420914774905315, +0.0789115157950706, +0.1479403360567813,
       +0.1479403360567813, +0.0789115157950706, +0.0789115157950706,
       +0.1479403360567813, +0.1479403360567813, +0.0789115157950706,
       +0.0420914774905315, +0.0789115157950706, +0.0789115157950706,
       +0.0420914774905315}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(QuadratureType::Cube, 1, 9, 5,
                                 {{{
                                      {{-0.9061798459386640}},
                                      {{-0.5384693101056831}},
                                      {{+0.0000000000000000}},
                                      {{+0.5384693101056831}},
                                      {{+0.9061798459386640}},
                                  }},
                                  {{+0.2369268850561891, +0.4786286704993665,
                                    +0.5688888888888889, +0.4786286704993665,
                                    +0.2369268850561891}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 2, 9, 25,
    {{{
         {{-0.9061798459386640, -0.9061798459386640}},
         {{-0.9061798459386640, -0.5384693101056831}},
         {{-0.9061798459386640, +0.0000000000000000}},
         {{-0.9061798459386640, +0.5384693101056831}},
         {{-0.9061798459386640, +0.9061798459386640}},
         {{-0.5384693101056831, -0.9061798459386640}},
         {{-0.5384693101056831, -0.5384693101056831}},
         {{-0.5384693101056831, +0.0000000000000000}},
         {{-0.5384693101056831, +0.5384693101056831}},
         {{-0.5384693101056831, +0.9061798459386640}},
         {{+0.0000000000000000, -0.9061798459386640}},
         {{+0.0000000000000000, -0.5384693101056831}},
         {{+0.0000000000000000, +0.0000000000000000}},
         {{+0.0000000000000000, +0.5384693101056831}},
         {{+0.0000000000000000, +0.9061798459386640}},
         {{+0.5384693101056831, -0.9061798459386640}},
         {{+0.5384693101056831, -0.5384693101056831}},
         {{+0.5384693101056831, +0.0000000000000000}},
         {{+0.5384693101056831, +0.5384693101056831}},
         {{+0.5384693101056831, +0.9061798459386640}},
         {{+0.9061798459386640, -0.9061798459386640}},
         {{+0.9061798459386640, -0.5384693101056831}},
         {{+0.9061798459386640, +0.0000000000000000}},
         {{+0.9061798459386640, +0.5384693101056831}},
         {{+0.9061798459386640, +0.9061798459386640}},
     }},
     {{+0.0561343488624286, +0.1134000000000000, +0.1347850723875209,
       +0.1134000000000000, +0.0561343488624286, +0.1134000000000000,
       +0.2290854042239911, +0.2722865325507507, +0.2290854042239911,
       +0.1134000000000000, +0.1347850723875209, +0.2722865325507507,
       +0.3236345679012346, +0.2722865325507507, +0.1347850723875209,
       +0.1134000000000000, +0.2290854042239911, +0.2722865325507507,
       +0.2290854042239911, +0.1134000000000000, +0.0561343488624286,
       +0.1134000000000000, +0.1347850723875209, +0.1134000000000000,
       +0.0561343488624286}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Cube, 3, 9, 125,
    {{{
         {{-0.9061798459386640, -0.9061798459386640, -0.9061798459386640}},
         {{-0.9061798459386640, -0.9061798459386640, -0.5384693101056831}},
         {{-0.9061798459386640, -0.9061798459386640, +0.0000000000000000}},
         {{-0.9061798459386640, -0.9061798459386640, +0.5384693101056831}},
         {{-0.9061798459386640, -0.9061798459386640, +0.9061798459386640}},
         {{-0.9061798459386640, -0.5384693101056831, -0.9061798459386640}},
         {{-0.9061798459386640, -0.5384693101056831, -0.5384693101056831}},
         {{-0.9061798459386640, -0.5384693101056831, +0.0000000000000000}},
         {{-0.9061798459386640, -0.5384693101056831, +0.5384693101056831}},
         {{-0.9061798459386640, -0.5384693101056831, +0.9061798459386640}},
         {{-0.9061798459386640, +0.0000000000000000, -0.9061798459386640}},
         {{-0.9061798459386640, +0.0000000000000000, -0.5384693101056831}},
         {{-0.9061798459386640, +0.0000000000000000, +0.0000000000000000}},
         {{-0.9061798459386640, +0.0000000000000000, +0.5384693101056831}},
         {{-0.9061798459386640, +0.0000000000000000, +0.9061798459386640}},
         {{-0.9061798459386640, +0.5384693101056831, -0.9061798459386640}},
         {{-0.9061798459386640, +0.5384693101056831, -0.5384693101056831}},
         {{-0.9061798459386640, +0.5384693101056831, +0.0000000000000000}},
         {{-0.9061798459386640, +0.5384693101056831, +0.5384693101056831}},
         {{-0.9061798459386640, +0.5384693101056831, +0.9061798459386640}},
         {{-0.9061798459386640, +0.9061798459386640, -0.9061798459386640}},
         {{-0.9061798459386640, +0.9061798459386640, -0.5384693101056831}},
         {{-0.9061798459386640, +0.9061798459386640, +0.0000000000000000}},
         {{-0.9061798459386640, +0.9061798459386640, +0.5384693101056831}},
         {{-0.9061798459386640, +0.9061798459386640, +0.9061798459386640}},
         {{-0.5384693101056831, -0.9061798459386640, -0.9061798459386640}},
         {{-0.5384693101056831, -0.9061798459386640, -0.5384693101056831}},
         {{-0.5384693101056831, -0.9061798459386640, +0.0000000000000000}},
         {{-0.5384693101056831, -0.9061798459386640, +0.5384693101056831}},
         {{-0.5384693101056831, -0.9061798459386640, +0.9061798459386640}},
         {{-0.5384693101056831, -0.5384693101056831, -0.9061798459386640}},
         {{-0.5384693101056831, -0.5384693101056831, -0.5384693101056831}},
         {{-0.5384693101056831, -0.5384693101056831, +0.0000000000000000}},
         {{-0.5384693101056831, -0.5384693101056831, +0.5384693101056831}},
         {{-0.5384693101056831, -0.5384693101056831, +0.9061798459386640}},
         {{-0.5384693101056831, +0.0000000000000000, -0.9061798459386640}},
         {{-0.5384693101056831, +0.0000000000000000, -0.5384693101056831}},
         {{-0.5384693101056831, +0.0000000000000000, +0.0000000000000000}},
         {{-0.5384693101056831, +0.0000000000000000, +0.5384693101056831}},
         {{-0.5384693101056831, +0.0000000000000000, +0.9061798459386640}},
         {{-0.5384693101056831, +0.5384693101056831, -0.9061798459386640}},
         {{-0.5384693101056831, +0.5384693101056831, -0.5384693101056831}},
         {{-0.5384693101056831, +0.5384693101056831, +0.0000000000000000}},
         {{-0.5384693101056831, +0.5384693101056831, +0.5384693101056831}},
         {{-0.5384693101056831, +0.5384693101056831, +0.9061798459386640}},
         {{-0.5384693101056831, +0.9061798459386640, -0.9061798459386640}},
         {{-0.5384693101056831, +0.9061798459386640, -0.5384693101056831}},
         {{-0.5384693101056831, +0.9061798459386640, +0.0000000000000000}},
         {{-0.5384693101056831, +0.9061798459386640, +0.5384693101056831}},
         {{-0.5384693101056831, +0.9061798459386640, +0.9061798459386640}},
         {{+0.0000000000000000, -0.9061798459386640, -0.9061798459386640}},
         {{+0.0000000000000000, -0.9061798459386640, -0.5384693101056831}},
         {{+0.0000000000000000, -0.9061798459386640, +0.0000000000000000}},
         {{+0.0000000000000000, -0.9061798459386640, +0.5384693101056831}},
         {{+0.0000000000000000, -0.9061798459386640, +0.9061798459386640}},
         {{+0.0000000000000000, -0.5384693101056831, -0.9061798459386640}},
         {{+0.0000000000000000, -0.5384693101056831, -0.5384693101056831}},
         {{+0.0000000000000000, -0.5384693101056831, +0.0000000000000000}},
         {{+0.0000000000000000, -0.5384693101056831, +0.5384693101056831}},
         {{+0.0000000000000000, -0.5384693101056831, +0.9061798459386640}},
         {{+0.0000000000000000, +0.0000000000000000, -0.9061798459386640}},
         {{+0.0000000000000000, +0.0000000000000000, -0.5384693101056831}},
         {{+0.0000000000000000, +0.0000000000000000, +0.0000000000000000}},
         {{+0.0000000000000000, +0.0000000000000000, +0.5384693101056831}},
         {{+0.0000000000000000, +0.0000000000000000, +0.9061798459386640}},
         {{+0.0000000000000000, +0.5384693101056831, -0.9061798459386640}},
         {{+0.0000000000000000, +0.5384693101056831, -0.5384693101056831}},
         {{+0.0000000000000000, +0.5384693101056831, +0.0000000000000000}},
         {{+0.0000000000000000, +0.5384693101056831, +0.5384693101056831}},
         {{+0.0000000000000000, +0.5384693101056831, +0.9061798459386640}},
         {{+0.0000000000000000, +0.9061798459386640, -0.9061798459386640}},
         {{+0.0000000000000000, +0.9061798459386640, -0.5384693101056831}},
         {{+0.0000000000000000, +0.9061798459386640, +0.0000000000000000}},
         {{+0.0000000000000000, +0.9061798459386640, +0.5384693101056831}},
         {{+0.0000000000000000, +0.9061798459386640, +0.9061798459386640}},
         {{+0.5384693101056831, -0.9061798459386640, -0.9061798459386640}},
         {{+0.5384693101056831, -0.9061798459386640, -0.5384693101056831}},
         {{+0.5384693101056831, -0.9061798459386640, +0.0000000000000000}},
         {{+0.5384693101056831, -0.9061798459386640, +0.5384693101056831}},
         {{+0.5384693101056831, -0.9061798459386640, +0.9061798459386640}},
         {{+0.5384693101056831, -0.5384693101056831, -0.9061798459386640}},
         {{+0.5384693101056831, -0.5384693101056831, -0.5384693101056831}},
         {{+0.5384693101056831, -0.5384693101056831, +0.0000000000000000}},
         {{+0.5384693101056831, -0.5384693101056831, +0.5384693101056831}},
         {{+0.5384693101056831, -0.5384693101056831, +0.9061798459386640}},
         {{+0.5384693101056831, +0.0000000000000000, -0.9061798459386640}},
         {{+0.5384693101056831, +0.0000000000000000, -0.5384693101056831}},
         {{+0.5384693101056831, +0.0000000000000000, +0.0000000000000000}},
         {{+0.5384693101056831, +0.0000000000000000, +0.5384693101056831}},
         {{+0.5384693101056831, +0.0000000000000000, +0.9061798459386640}},
         {{+0.5384693101056831, +0.5384693101056831, -0.9061798459386640}},
         {{+0.5384693101056831, +0.5384693101056831, -0.5384693101056831}},
         {{+0.5384693101056831, +0.5384693101056831, +0.0000000000000000}},
         {{+0.5384693101056831, +0.5384693101056831, +0.5384693101056831}},
         {{+0.5384693101056831, +0.5384693101056831, +0.9061798459386640}},
         {{+0.5384693101056831, +0.9061798459386640, -0.9061798459386640}},
         {{+0.5384693101056831, +0.9061798459386640, -0.5384693101056831}},
         {{+0.5384693101056831, +0.9061798459386640, +0.0000000000000000}},
         {{+0.5384693101056831, +0.9061798459386640, +0.5384693101056831}},
         {{+0.5384693101056831, +0.9061798459386640, +0.9061798459386640}},
         {{+0.9061798459386640, -0.9061798459386640, -0.9061798459386640}},
         {{+0.9061798459386640, -0.9061798459386640, -0.5384693101056831}},
         {{+0.9061798459386640, -0.9061798459386640, +0.0000000000000000}},
         {{+0.9061798459386640, -0.9061798459386640, +0.5384693101056831}},
         {{+0.9061798459386640, -0.9061798459386640, +0.9061798459386640}},
         {{+0.9061798459386640, -0.5384693101056831, -0.9061798459386640}},
         {{+0.9061798459386640, -0.5384693101056831, -0.5384693101056831}},
         {{+0.9061798459386640, -0.5384693101056831, +0.0000000000000000}},
         {{+0.9061798459386640, -0.5384693101056831, +0.5384693101056831}},
         {{+0.9061798459386640, -0.5384693101056831, +0.9061798459386640}},
         {{+0.9061798459386640, +0.0000000000000000, -0.9061798459386640}},
         {{+0.9061798459386640, +0.0000000000000000, -0.5384693101056831}},
         {{+0.9061798459386640, +0.0000000000000000, +0.0000000000000000}},
         {{+0.9061798459386640, +0.0000000000000000, +0.5384693101056831}},
         {{+0.9061798459386640, +0.0000000000000000, +0.9061798459386640}},
         {{+0.9061798459386640, +0.5384693101056831, -0.9061798459386640}},
         {{+0.9061798459386640, +0.5384693101056831, -0.5384693101056831}},
         {{+0.9061798459386640, +0.5384693101056831, +0.0000000000000000}},
         {{+0.9061798459386640, +0.5384693101056831, +0.5384693101056831}},
         {{+0.9061798459386640, +0.5384693101056831, +0.9061798459386640}},
         {{+0.9061798459386640, +0.9061798459386640, -0.9061798459386640}},
         {{+0.9061798459386640, +0.9061798459386640, -0.5384693101056831}},
         {{+0.9061798459386640, +0.9061798459386640, +0.0000000000000000}},
         {{+0.9061798459386640, +0.9061798459386640, +0.5384693101056831}},
         {{+0.9061798459386640, +0.9061798459386640, +0.9061798459386640}},
     }},
     {{+0.0132997364206326, +0.0268675087653718, +0.0319342073528483,
       +0.0268675087653718, +0.0132997364206326, +0.0268675087653718,
       +0.0542764912346282, +0.0645120000000000, +0.0542764912346282,
       +0.0268675087653718, +0.0319342073528483, +0.0645120000000000,
       +0.0766777300693452, +0.0645120000000000, +0.0319342073528483,
       +0.0268675087653718, +0.0542764912346282, +0.0645120000000000,
       +0.0542764912346282, +0.0268675087653718, +0.0132997364206326,
       +0.0268675087653718, +0.0319342073528483, +0.0268675087653718,
       +0.0132997364206326, +0.0268675087653718, +0.0542764912346282,
       +0.0645120000000000, +0.0542764912346282, +0.0268675087653718,
       +0.0542764912346282, +0.1096468424545388, +0.1303241410696483,
       +0.1096468424545388, +0.0542764912346282, +0.0645120000000000,
       +0.1303241410696483, +0.1549007829622049, +0.1303241410696483,
       +0.0645120000000000, +0.0542764912346282, +0.1096468424545388,
       +0.1303241410696483, +0.1096468424545388, +0.0542764912346282,
       +0.0268675087653718, +0.0542764912346282, +0.0645120000000000,
       +0.0542764912346282, +0.0268675087653718, +0.0319342073528483,
       +0.0645120000000000, +0.0766777300693452, +0.0645120000000000,
       +0.0319342073528483, +0.0645120000000000, +0.1303241410696483,
       +0.1549007829622049, +0.1303241410696483, +0.0645120000000000,
       +0.0766777300693452, +0.1549007829622049, +0.1841121097393690,
       +0.1549007829622049, +0.0766777300693452, +0.0645120000000000,
       +0.1303241410696483, +0.1549007829622049, +0.1303241410696483,
       +0.0645120000000000, +0.0319342073528483, +0.0645120000000000,
       +0.0766777300693452, +0.0645120000000000, +0.0319342073528483,
       +0.0268675087653718, +0.0542764912346282, +0.0645120000000000,
       +0.0542764912346282, +0.0268675087653718, +0.0542764912346282,
       +0.1096468424545388, +0.1303241410696483, +0.1096468424545388,
       +0.0542764912346282, +0.0645120000000000, +0.1303241410696483,
       +0.1549007829622049, +0.1303241410696483, +0.0645120000000000,
       +0.0542764912346282, +0.1096468424545388, +0.1303241410696483,
       +0.1096468424545388, +0.0542764912346282, +0.0268675087653718,
       +0.0542764912346282, +0.0645120000000000, +0.0542764912346282,
       +0.0268675087653718, +0.0132997364206326, +0.0268675087653718,
       +0.0319342073528483, +0.0268675087653718, +0.0132997364206326,
       +0.0268675087653718, +0.0542764912346282, +0.0645120000000000,
       +0.0542764912346282, +0.0268675087653718, +0.0319342073528483,
       +0.0645120000000000, +0.0766777300693452, +0.0645120000000000,
       +0.0319342073528483, +0.0268675087653718, +0.0542764912346282,
       +0.0645120000000000, +0.0542764912346282, +0.0268675087653718,
       +0.0132997364206326, +0.0268675087653718, +0.0319342073528483,
       +0.0268675087653718, +0.0132997364206326}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Simplex, 2, 1, 1,
    {{{
         {{+0.3333333333333333, +0.3333333333333333}},
     }},
     {{+0.5000000000000000}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Simplex, 2, 2, 3,
    {{{
         {{+0.6666666666666667, +0.1666666666666667}},
         {{+0.1666666666666667, +0.1666666666666667}},
         {{+0.1666666666666667, +0.6666666666666667}},
     }},
     {{+0.1666666666666667, +0.1666666666666667, +0.1666666666666667}}});

AE108_ELEMENTS_QUADRATURE_DEFINE(
    QuadratureType::Simplex, 2, 3, 4,
    {{{
         {{+0.3333333333333333, +0.3333333333333333}},
         {{+0.6000000000000000, +0.2000000000000000}},
         {{+0.2000000000000000, +0.2000000000000000}},
         {{+0.2000000000000000, +0.6000000000000000}},
     }},
     {{-0.2812500000000000, +0.2604166666666667, +0.2604166666666667,
       +0.2604166666666667}}});

template <QuadratureType Type_, std::size_t Dimension_, std::size_t Order_>
struct IntegrateTrait<Quadrature<Type_, Dimension_, Order_>> {
  template <class Quadrature, class R, class F, class... Args>
  typename std::decay<R>::type operator()(
      F &&f, R &&init,
      const typename Quadrature::template Collection<Args> &... args) const
      noexcept {
    return eval(std::forward<F>(f), std::forward<R>(init),
                typename Quadrature::size_type{0},
                Quadrature::data.points.begin(), Quadrature::data.points.end(),
                Quadrature::data.weights.begin(), args.begin()...);
  }

private:
  template <class R, class F, class ID, class PointIterator,
            class WeightIterator, class... Iterators>
  static typename std::decay<R>::type
  eval(F f, R init, ID id, PointIterator point, const PointIterator point_end,
       WeightIterator weight, Iterators... iterators) noexcept {
    while (point != point_end) {
      init += *(weight++) * f(id++, *(point++), (*(iterators++))...);
    }
    return init;
  }
};

} // namespace quadrature
} // namespace elements
} // namespace ae108