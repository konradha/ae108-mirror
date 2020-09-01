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
void resizeIfPossible(std::array<T, Size> *const container,
                      const typename std::array<T, Size>::size_type size) {
  assert(container);
  assert(Size == size);
}
} // namespace utilities
} // namespace assembly
} // namespace ae108