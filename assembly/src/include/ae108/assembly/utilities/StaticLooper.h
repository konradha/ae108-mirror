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

#include <cstddef>
#include <utility>

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief Executes functor with parameters value.get<N>() and args for all N
 * from From to To (excluding).
 */
template <std::size_t To, std::size_t From = 0u> struct StaticLooper {
  template <class Functor, class Value, class... Args>
  static void run(Functor &&functor, Value &&value, Args &&...args);
};

/**
 * @brief Special case: Does not execute functor.
 */
template <std::size_t To> struct StaticLooper<To, To> {
  template <class Functor, class Value, class... Args>
  static void run(Functor &&, Value &&, Args &&...args);
};
} // namespace utilities
} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace assembly {
namespace utilities {

template <std::size_t To, std::size_t From>
template <class Functor, class Value, class... Args>
void StaticLooper<To, From>::run(Functor &&functor, Value &&value,
                                 Args &&...args) {
  functor(value.template get<From>(), args...);
  StaticLooper<To, From + 1>::run(std::forward<Functor>(functor),
                                  std::forward<Value>(value),
                                  std::forward<Args>(args)...);
}

template <std::size_t To>
template <class Functor, class Value, class... Args>
void StaticLooper<To, To>::run(Functor &&, Value &&, Args &&...) {}
} // namespace utilities
} // namespace assembly
} // namespace ae108