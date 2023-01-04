// © 2021 ETH Zurich, Mechanics and Materials Lab
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

// see G. Dhatt, G. Touzot, E. Lefrançois, "Finite Element Method"
// ISTE Ltd, London, 2012, p. 139.
struct Tet4 : ShapeBase<std::size_t, double, 3, 4> {};

AE108_ELEMENTS_SHAPE_DEFINE_VALUE(Tet4, {{
                                            1 - xi[0] - xi[1] - xi[2],
                                            xi[0],
                                            xi[1],
                                            xi[2],
                                        }});

AE108_ELEMENTS_SHAPE_DEFINE_GRADIENTS(Tet4, {{
                                                {{-1., -1., -1.}},
                                                {{1., 0., 0.}},
                                                {{0., 1., 0.}},
                                                {{0., 0., 1.}},
                                            }});

AE108_ELEMENTS_SHAPE_DEFINE_POINTS(Tet4, {{
                                             {{0., 0., 0.}},
                                             {{1., 0., 0.}},
                                             {{0., 1., 0.}},
                                             {{0., 0., 1.}},
                                         }});

} // namespace shape
} // namespace elements
} // namespace ae108