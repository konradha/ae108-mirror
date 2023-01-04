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

#include <boost/version.hpp>

#if BOOST_VERSION >= 107000
#include <boost/math/differentiation/finite_difference.hpp>
#else
#include <boost/math/tools/numerical_differentiation.hpp>
#endif
#include <type_traits>

namespace ae108 {
namespace elements {
namespace tensor {

/**
 * @brief Calls a Boost library function to perform finite difference
 * differentiation. Only supports differentiation of real-valued functions of
 * one real-valued parameter.
 */
template <class T, class F>
typename std::decay<T>::type differentiate(F &&f, T &&t) noexcept {
  return
#if BOOST_VERSION >= 107000
      boost::math::differentiation::finite_difference_derivative
#else
      boost::math::tools::finite_difference_derivative
#endif
      (std::forward<F>(f), std::forward<T>(t));
}

} // namespace tensor
} // namespace elements
} // namespace ae108
