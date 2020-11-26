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
                   const Viewer<Policy> &viewer);

/**
 * @brief Writes a vector to a Viewer. The concrete format depends on the file
 * extension (vtk/vtu).
 */
template <class Policy>
void writeToViewer(const distributed<Vector<Policy>> &data,
                   const distributed<Vector<Policy>> &coordinates,
                   const Viewer<Policy> &viewer);

extern template void writeToViewer<SequentialComputePolicy>(
    const Mesh<SequentialComputePolicy> &,
    const distributed<Vector<SequentialComputePolicy>> &,
    const Viewer<SequentialComputePolicy> &);
extern template void writeToViewer<ParallelComputePolicy>(
    const Mesh<ParallelComputePolicy> &,
    const distributed<Vector<ParallelComputePolicy>> &,
    const Viewer<ParallelComputePolicy> &);
extern template void writeToViewer<SequentialComputePolicy>(
    const distributed<Vector<SequentialComputePolicy>> &,
    const distributed<Vector<SequentialComputePolicy>> &,
    const Viewer<SequentialComputePolicy> &);
extern template void writeToViewer<ParallelComputePolicy>(
    const distributed<Vector<ParallelComputePolicy>> &,
    const distributed<Vector<ParallelComputePolicy>> &,
    const Viewer<ParallelComputePolicy> &);

} // namespace cpppetsc
} // namespace ae108

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

template <class Policy>
void writeToViewer(const Mesh<Policy> &mesh,
                   const distributed<Vector<Policy>> &coordinates,
                   const Viewer<Policy> &viewer) {
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

  Policy::handleError(DMView(dm.get(), viewer.data()));
}

template <class Policy>
void writeToViewer(const distributed<Vector<Policy>> &data,
                   const distributed<Vector<Policy>> &coordinates,
                   const Viewer<Policy> &viewer) {
  const auto getDM = [](const Vector<Policy> &x) {
    auto dm = DM{};
    Policy::handleError(VecGetDM(x.data(), &dm));
    assert(dm);
    return UniqueEntity<DM>(dm, [](DM) {});
  };

  // clone DM that defines the layout of the data vector
  // to be able to set the coordinate DM
  const auto dm = [&]() {
    auto dm = DM{};
    const auto reference = getDM(data).get();
    Policy::handleError(DMClone(reference, &dm));
    return makeUniqueEntity<Policy>(dm);
  }();

  // configure the default section to be able to create
  // a copy of the data vector
  auto section = PetscSection{};
  Policy::handleError(DMGetDefaultSection(getDM(data).get(), &section));
  Policy::handleError(DMSetDefaultSection(dm.get(), section));

  // set the coordinate DM which is required by VecView
  const auto coordinateDM = getDM(coordinates);
  Policy::handleError(DMSetCoordinateDM(dm.get(), coordinateDM.get()));
  Policy::handleError(DMSetCoordinates(dm.get(), coordinates.unwrap().data()));

  // copy the data to the vector with the correct DM
  const auto vec = [&]() {
    auto vec = Vec{};
    Policy::handleError(DMCreateGlobalVector(dm.get(), &vec));
    return makeUniqueEntity<Policy>(vec);
  }();
  Policy::handleError(VecCopy(data.unwrap().data(), vec.get()));

  Policy::handleError(VecView(vec.get(), viewer.data()));
}

} // namespace cpppetsc
} // namespace ae108
