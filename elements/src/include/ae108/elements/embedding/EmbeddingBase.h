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

#include <cstddef>

namespace ae108 {
namespace elements {
namespace embedding {

template <std::size_t FromDimension_, class ReferencePoint_,
          std::size_t ToDimension_, class PhysicalPoint_,
          class SizeType_ = std::size_t, class ValueType_ = double>
struct EmbeddingBase {
  using value_type = ValueType_;
  using size_type = SizeType_;

  /**
   * @brief The dimension of reference space.
   */
  static constexpr size_type reference_dimension() noexcept {
    return FromDimension_;
  }

  /**
   * @brief The dimension of physical space.
   */
  static constexpr size_type physical_dimension() noexcept {
    return ToDimension_;
  }

  using ReferencePoint = ReferencePoint_;
  using PhysicalPoint = PhysicalPoint_;

protected:
  ~EmbeddingBase() = default;
};

} // namespace embedding
} // namespace elements
} // namespace ae108