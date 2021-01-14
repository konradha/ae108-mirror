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

#include <array>
#include <cstddef>

namespace ae108 {
namespace elements {
namespace tensor {

namespace detail {
/**
 * @brief A helper class to create an array of arrays from a list of sizes.
 */
template <class, std::size_t...> struct TensorImpl;

template <class T_, std::size_t Size_> struct TensorImpl<T_, Size_> {
  using type = std::array<T_, Size_>;
};

template <class T_, std::size_t Size_, std::size_t... Sizes_>
struct TensorImpl<T_, Size_, Sizes_...> {
  using type = std::array<typename TensorImpl<T_, Sizes_...>::type, Size_>;
};
} // namespace detail

/**
 * @brief An alias for an array of arrays. For instance, in the case of a
 * matrix, the first size parameter is the number of rows and the second size
 * parameter is the number of columns.
 *
 * @tparam T The value type.
 *
 * @remark The values are stored in row-major format. For instance, a matrix is
 * stored as an array of rows.
 */
template <class T, std::size_t... Sizes>
using Tensor = typename detail::TensorImpl<T, Sizes...>::type;

} // namespace tensor
} // namespace elements
} // namespace ae108