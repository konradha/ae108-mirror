// © 2022 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/integrator/IntegrateShapeTrait.h"
#include <type_traits>
#include <utility>

namespace ae108 {
namespace elements {
namespace integrator {

/**
 * @brief Integrates the function numerically. Passes 2 parameters at every
 * quadrature point: the index of the quadrature point, the value of the shape
 * functions at the quadrature point.
 */
template <class Integrator, class R, class F>
std::decay_t<R> integrate_shape(const Integrator &integrator, F &&f,
                                R &&init) noexcept {
  return IntegrateShapeTrait<Integrator>().template operator()<Integrator>(
      integrator, std::forward<F>(f), std::forward<R>(init));
}

} // namespace integrator
} // namespace elements
} // namespace ae108