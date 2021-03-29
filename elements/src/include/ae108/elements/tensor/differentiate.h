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
