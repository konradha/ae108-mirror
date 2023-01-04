// © 2022 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/integrator/ComputeVolumeTrait.h"

namespace ae108 {
namespace elements {
namespace integrator {

/**
 * @brief Integrates the constant function 1.
 */
template <class Integrator>
typename Integrator::real_type volume(const Integrator &integrator) noexcept {
  return ComputeVolumeTrait<Integrator>().template operator()(integrator);
}

} // namespace integrator
} // namespace elements
} // namespace ae108