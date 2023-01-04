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

#include "ae108/elements/integrator/IntegrateTrait.h"
#include <type_traits>
#include <utility>

namespace ae108 {
namespace elements {
namespace integrator {

/**
 * @tparam DiscretizedFunction Equal to an Integrator::DiscretizedFunction<...>.
 */
template <class Integrator, class R, class F, class DiscretizedFunction>
typename std::decay<R>::type integrate(const Integrator &integrator, F &&f,
                                       const DiscretizedFunction &u,
                                       R &&init) noexcept {
  return IntegrateTrait<Integrator>().template operator()<Integrator>(
      integrator, std::forward<F>(f), u, std::forward<R>(init));
}

} // namespace integrator
} // namespace elements
} // namespace ae108