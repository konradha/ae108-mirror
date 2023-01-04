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
          const typename Quadrature::template Collection<Args> &...args) {
  return IntegrateTrait<Quadrature>().template operator()<Quadrature>(
      std::forward<F>(f), std::forward<R>(init), args...);
}

} // namespace quadrature
} // namespace elements
} // namespace ae108