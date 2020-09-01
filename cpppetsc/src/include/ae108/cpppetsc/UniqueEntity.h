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

#include <functional>
#include <memory>
#include <petscsys.h>
#include <type_traits>

namespace ae108 {
namespace cpppetsc {

template <typename T>
using UniqueEntity = std::unique_ptr<typename std::remove_pointer<T>::type,
                                     std::function<void(T)>>;

namespace detail {
template <class T> PetscErrorCode callDestructor(T *const ptr) noexcept;
}

/**
 * @brief Create a UniqueEntity with the default destructor.
 *
 * @tparam Policy The policy used to handle errors.
 */
template <typename Policy, typename T>
UniqueEntity<T> makeUniqueEntity(const T &t) {
  return UniqueEntity<T>{
      t, [](T t) { Policy::handleError(detail::callDestructor(&t)); }};
}
} // namespace cpppetsc
} // namespace ae108