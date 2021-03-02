// © 2021 ETH Zurich, Mechanics and Materials Lab
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

// see G. Dhatt, G. Touzot, E. Lefrançois, "Finite Element Method"
// ISTE Ltd, London, 2012, p. 140.
struct Tet10 : ShapeBase<std::size_t, double, 3, 10> {};

AE108_ELEMENTS_SHAPE_DEFINE_VALUE(
    Tet10,
    {{
        -(1. - xi[0] - xi[1] - xi[2]) * (1. - 2. * (1 - xi[0] - xi[1] - xi[2])),
        -xi[0] * (1 - 2. * xi[0]),
        -xi[1] * (1. - 2. * xi[1]),
        -xi[2] * (1. - 2. * xi[2]),
        4. * xi[0] * (1 - xi[0] - xi[1] - xi[2]),
        4. * xi[0] * xi[1],
        4. * xi[1] * (1 - xi[0] - xi[1] - xi[2]),
        4. * xi[2] * (1 - xi[0] - xi[1] - xi[2]),
        4. * xi[0] * xi[2],
        4. * xi[1] * xi[2],
    }});

AE108_ELEMENTS_SHAPE_DEFINE_GRADIENTS(
    Tet10,
    {{
        {{1. - 4. * (1 - xi[0] - xi[1] - xi[2]),
          1. - 4. * (1 - xi[0] - xi[1] - xi[2]),
          1. - 4. * (1 - xi[0] - xi[1] - xi[2])}},
        {{-1. + 4. * xi[0], 0., 0.}},
        {{0., -1. + 4. * xi[1], 0.}},
        {{0., 0., -1. + 4. * xi[2]}},
        {{4. * ((1 - xi[0] - xi[1] - xi[2]) - xi[0]), -4. * xi[0],
          -4. * xi[0]}},
        {{4. * xi[1], 4. * xi[0], 0.}},
        {{-4. * xi[1], 4. * ((1 - xi[0] - xi[1] - xi[2]) - xi[1]), -4 * xi[1]}},
        {{-4. * xi[2], -4. * xi[2],
          4. * ((1 - xi[0] - xi[1] - xi[2]) - xi[2])}},
        {{4. * xi[2], 0., 4. * xi[0]}},
        {{0., 4. * xi[2], 4. * xi[1]}},
    }});

AE108_ELEMENTS_SHAPE_DEFINE_POINTS(Tet10, {{
                                              {{0., 0., 0.}},
                                              {{1., 0., 0.}},
                                              {{0., 1., 0.}},
                                              {{0., 0., 1.}},
                                              {{.5, 0., 0.}},
                                              {{.5, .5, 0.}},
                                              {{0., .5, 0.}},
                                              {{0., 0., .5}},
                                              {{.5, 0., .5}},
                                              {{0., .5, .5}},
                                          }});

} // namespace shape
} // namespace elements
} // namespace ae108