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

#include "ae108/elements/integrator/VolumeTrait.h"

namespace ae108 {
namespace elements {
namespace integrator {

/**
 * @brief Integrates the constant function 1.
 */
template <class Integrator>
typename Integrator::real_type volume(const Integrator &integrator) noexcept {
  return VolumeTrait<Integrator>().template operator()(integrator);
}

} // namespace integrator
} // namespace elements
} // namespace ae108