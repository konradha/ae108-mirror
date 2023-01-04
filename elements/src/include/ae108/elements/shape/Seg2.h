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

struct Seg2 : ShapeBase<std::size_t, double, 1, 2> {};

AE108_ELEMENTS_SHAPE_DEFINE_POINTS(Seg2, {{
                                             {{-1.}},
                                             {{1.}},
                                         }});

AE108_ELEMENTS_SHAPE_DEFINE_VALUE(Seg2, {{
                                            -.5 * (xi[0] - 1.),
                                            .5 * (xi[0] + 1.),
                                        }});

AE108_ELEMENTS_SHAPE_DEFINE_GRADIENTS(Seg2, {{
                                                {{-.5}},
                                                {{.5}},
                                            }});

} // namespace shape
} // namespace elements
} // namespace ae108