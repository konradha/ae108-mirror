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

#include "ae108/elements/quadrature/Quadrature.h"

#define AE108_ELEMENTS_QUADRATURE_DEFINE_CC(type, dimension, order)            \
  constexpr typename Quadrature<type, dimension, order>::Data                  \
      Quadrature<type, dimension, order>::data

namespace ae108 {
namespace elements {
namespace quadrature {

AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 3);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 3);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 3);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 5);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 5);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 5);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 7);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 7);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 7);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 1, 9);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 2, 9);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Cube, 3, 9);

AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Simplex, 2, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Simplex, 2, 2);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Simplex, 2, 3);

AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Simplex, 3, 1);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Simplex, 3, 2);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Simplex, 3, 3);
AE108_ELEMENTS_QUADRATURE_DEFINE_CC(QuadratureType::Simplex, 3, 4);

} // namespace quadrature
} // namespace elements
} // namespace ae108