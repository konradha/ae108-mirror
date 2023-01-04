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

#include "ae108/cpppetsc/TaggedEntity.h"

namespace ae108 {
namespace cpppetsc {

struct LocalTag;
struct DistributedTag;
struct GlobalTag;

/**
 * @brief A tag to annotate a "local" vector. This vector only contains the data
 * for all local vertices.
 */
template <class Entity> using local = TaggedEntity<Entity, LocalTag>;

/**
 * @brief A tag to annotate a "distributed" vector. This vector contains the
 * data for all vertices, but the data is distributed between MPI nodes.
 */
template <class Entity>
using distributed = TaggedEntity<Entity, DistributedTag>;

/**
 * @brief A tag to annotate a "global" vector. This vector contains the
 * data for all vertices, and all data is available locally.
 */
template <class Entity> using global = TaggedEntity<Entity, GlobalTag>;

} // namespace cpppetsc
} // namespace ae108