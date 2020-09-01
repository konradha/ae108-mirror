// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

namespace ae108 {
namespace solve {
namespace dynamics {

template <class VectorType> struct DynamicState {
  using value_type = VectorType;

  value_type displacements;
  value_type velocities;
  value_type accelerations;
};

} // namespace dynamics
} // namespace solve
} // namespace ae108
