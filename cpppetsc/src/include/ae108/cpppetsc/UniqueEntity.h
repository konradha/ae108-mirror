// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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