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

#define AE108_ELEMENTS_QUADRATURE_DEFINE_CC(type, dimension, order)            \
  constexpr typename Quadrature<type, dimension, order>::Data                  \
      Quadrature<type, dimension, order>::data;

namespace ae108 {
namespace elements {
namespace quadrature {

AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 2);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 2);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 2);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 3);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 3);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 3);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 4);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 4);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 4);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 5);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 5);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 5);

} // namespace quadrature
} // namespace elements
} // namespace ae108