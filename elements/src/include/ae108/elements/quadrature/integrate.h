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
#include <utility>

namespace ae108 {
namespace elements {
namespace quadrature {

/**
 * @brief Integrates the function f and adds the result to init. The function f
 * is called with the parameters id, point, args_1, args_2, ... . Here, "id" is
 * the quadrature point id, "point" is the quadrature point, and "args_1" is the
 * element of the first args corresponding to the quadrature point id.
 */
template <class Quadrature, class R, class F, class... Args>
constexpr typename std::decay<R>::type
integrate(F &&f, R &&init,
          const typename Quadrature::template Collection<Args> &... args) {
  return IntegrateTrait<Quadrature>().template operator()<Quadrature>(
      std::forward<F>(f), std::forward<R>(init), args...);
}

} // namespace quadrature
} // namespace elements
} // namespace ae108