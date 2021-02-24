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

#include "ae108/elements/shape/GradientTrait.h"
#include "ae108/elements/shape/PointTrait.h"
#include "ae108/elements/shape/ShapeBase.h"
#include "ae108/elements/shape/ValueTrait.h"
#include <cstddef>

namespace ae108 {
namespace elements {
namespace shape {

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