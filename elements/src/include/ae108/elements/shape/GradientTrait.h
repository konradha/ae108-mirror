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

namespace ae108 {
namespace elements {
namespace shape {

template <class Shape> struct GradientTrait;

#define AE108_ELEMENTS_SHAPE_DEFINE_GRADIENTS(name, ...)                       \
  template <> struct GradientTrait<name> {                                     \
    constexpr typename name::template Collection<typename name::Gradient>      \
    operator()(const name &, const typename name::Point &xi) noexcept {        \
      return __VA_ARGS__;                                                      \
    }                                                                          \
  }

} // namespace shape
} // namespace elements
} // namespace ae108