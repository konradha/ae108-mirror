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

#include "ae108/elements/embedding/AutomaticJacobianTrait.h"

namespace ae108 {
namespace elements {
namespace embedding {

/**
 * @brief Computes the Jacobian using numerical differentiation of the
 * embedding.
 */
template <class Embedding>
typename AutomaticJacobianTrait<Embedding>::template Jacobian<Embedding>
automatic_jacobian(const Embedding &embedding,
                   const typename Embedding::ReferencePoint &point) {
  return AutomaticJacobianTrait<Embedding>()(embedding, point);
}

} // namespace embedding
} // namespace elements
} // namespace ae108
