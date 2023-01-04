// © 2020 ETH Zurich, Mechanics and Materials Lab
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