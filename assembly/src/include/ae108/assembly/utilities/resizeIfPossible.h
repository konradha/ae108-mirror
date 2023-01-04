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

#include <array>
#include <cstddef>
#include <vector>

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief Resizes the vector to size.
 * @param container Valid nonzero pointer.
 */
template <class T>
void resizeIfPossible(std::vector<T> *const container,
                      const typename std::vector<T>::size_type size);

/**
 * @brief Asserts that the array has the correct size.
 * @param container Valid nonzero pointer.
 */
template <class T, std::size_t Size>
void resizeIfPossible(std::array<T, Size> *const container,
                      const typename std::array<T, Size>::size_type size);
} // namespace utilities
} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include <cassert>

namespace ae108 {
namespace assembly {
namespace utilities {

template <class T>
void resizeIfPossible(std::vector<T> *const container,
                      const typename std::vector<T>::size_type size) {
  assert(container);
  container->resize(size);
}

template <class T, std::size_t Size>
void resizeIfPossible(std::array<T, Size> *const container [[maybe_unused]],
                      const typename std::array<T, Size>::size_type size
                      [[maybe_unused]]) {
  assert(container);
  assert(Size == size);
}
} // namespace utilities
} // namespace assembly
} // namespace ae108