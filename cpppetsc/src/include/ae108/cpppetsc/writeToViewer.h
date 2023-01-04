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

#include "ae108/cpppetsc/Mesh_fwd.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector_fwd.h"
#include "ae108/cpppetsc/Viewer.h"
#include <cassert>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Writes a mesh to a Viewer.
 */
template <class Policy>
void writeToViewer(const Mesh<Policy> &mesh,
                   const distributed<Vector<Policy>> &coordinates,
                   Viewer<Policy> *viewer);

/**
 * @brief Writes a vector to a Viewer.
 */
template <class Policy>
void writeToViewer(const distributed<Vector<Policy>> &data,
                   Viewer<Policy> *viewer);

extern template void writeToViewer<SequentialComputePolicy>(
    const Mesh<SequentialComputePolicy> &,
    const distributed<Vector<SequentialComputePolicy>> &,
    Viewer<SequentialComputePolicy> *);
extern template void writeToViewer<ParallelComputePolicy>(
    const Mesh<ParallelComputePolicy> &,
    const distributed<Vector<ParallelComputePolicy>> &,
    Viewer<ParallelComputePolicy> *);
extern template void writeToViewer<SequentialComputePolicy>(
    const distributed<Vector<SequentialComputePolicy>> &,
    Viewer<SequentialComputePolicy> *);
extern template void writeToViewer<ParallelComputePolicy>(
    const distributed<Vector<ParallelComputePolicy>> &,
    Viewer<ParallelComputePolicy> *);

} // namespace cpppetsc
} // namespace ae108

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

template <class Policy>
void writeToViewer(const Mesh<Policy> &mesh,
                   const distributed<Vector<Policy>> &coordinates,
                   Viewer<Policy> *const viewer) {
  assert(viewer);

  const auto coordinateDM = [&]() {
    auto dm = DM{};
    Policy::handleError(VecGetDM(coordinates.unwrap().data(), &dm));
    assert(dm);
    return UniqueEntity<DM>(dm, [](DM) {});
  }();

  // copy the DM to be able to set the coordinate DM, which is
  // required by DMView
  auto dm = [&]() {
    auto dm = DM{};
    Policy::handleError(DMClone(mesh.data(), &dm));
    return makeUniqueEntity<Policy>(dm);
  }();
  Policy::handleError(DMSetCoordinateDM(dm.get(), coordinateDM.get()));
  Policy::handleError(DMSetCoordinates(dm.get(), coordinates.unwrap().data()));

  Policy::handleError(DMView(dm.get(), viewer->data()));
}

template <class Policy>
void writeToViewer(const distributed<Vector<Policy>> &data,
                   Viewer<Policy> *const viewer) {
  assert(viewer);
  Policy::handleError(VecView(data.unwrap().data(), viewer->data()));
}

} // namespace cpppetsc
} // namespace ae108
