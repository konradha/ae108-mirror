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

#pragma once

#include "ae108/elements/shape/GradientTrait.h"
#include "ae108/elements/shape/PointTrait.h"
#include "ae108/elements/shape/ShapeBase.h"
#include "ae108/elements/shape/ValueTrait.h"
#include <cstddef>

namespace ae108 {
namespace elements {
namespace shape {

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.265
struct Tri6 : ShapeBase<std::size_t, double, 2, 6> {};

AE108_ELEMENTS_SHAPE_DEFINE_VALUE(Tri6, {{
                                            (1 - xi[0] - xi[1]) *
                                                (1 - 2 * xi[0] - 2 * xi[1]),
                                            xi[0] * (2 * xi[0] - 1),
                                            xi[1] * (2 * xi[1] - 1),
                                            4 * xi[0] * (1 - xi[0] - xi[1]),
                                            4 * xi[0] * xi[1],
                                            4 * xi[1] * (1 - xi[0] - xi[1]),
                                        }});

AE108_ELEMENTS_SHAPE_DEFINE_GRADIENTS(
    Tri6, {{
              {{4. * (xi[0] + xi[1]) - 3., 4. * (xi[0] + xi[1]) - 3.}},
              {{4. * xi[0] - 1., 0.}},
              {{0., 4 * xi[1] - 1}},
              {{4. * (1 - 2 * xi[0] - xi[1]), -4. * xi[0]}},
              {{4. * xi[1], 4. * xi[0]}},
              {{-4. * xi[1], 4. * (1 - xi[0] - 2. * xi[1])}},
          }});

AE108_ELEMENTS_SHAPE_DEFINE_POINTS(Tri6, {{
                                             {{0., 0.}},
                                             {{1., 0.}},
                                             {{0., 1.}},
                                             {{0.5, 0.}},
                                             {{0.5, 0.5}},
                                             {{0., 0.5}},
                                         }});

} // namespace shape
} // namespace elements
} // namespace ae108