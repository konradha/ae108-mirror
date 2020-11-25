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

#include "ae108/cpppetsc/InvalidParametersException.h"
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
 *
 * @throw InvalidParametersException if the coordinates of `mesh` are not set.
 */
template <class Policy>
void writeToViewer(const Mesh<Policy> &mesh, const Viewer<Policy> &viewer);

/**
 * @brief Writes a vector to a Viewer.
 *
 * @throw InvalidParametersException if the `data` vector has no associated mesh
 * or if the coordinates of this mesh are not set.
 */
template <class Policy>
void writeToViewer(const distributed<Vector<Policy>> &data,
                   const Viewer<Policy> &viewer);

extern template void
writeToViewer<SequentialComputePolicy>(const Mesh<SequentialComputePolicy> &,
                                       const Viewer<SequentialComputePolicy> &);
extern template void
writeToViewer<ParallelComputePolicy>(const Mesh<ParallelComputePolicy> &,
                                     const Viewer<ParallelComputePolicy> &);
extern template void writeToViewer<SequentialComputePolicy>(
    const distributed<Vector<SequentialComputePolicy>> &,
    const Viewer<SequentialComputePolicy> &);
extern template void writeToViewer<ParallelComputePolicy>(
    const distributed<Vector<ParallelComputePolicy>> &,
    const Viewer<ParallelComputePolicy> &);

} // namespace cpppetsc
} // namespace ae108

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {
namespace detail {
template <class Policy> void throwIfNoCoordinateDM(const DM dm) {
  auto coordinateDM = DM{};
  Policy::handleError(DMGetCoordinateDM(dm, &coordinateDM));
  if (!coordinateDM) {
    throw InvalidParametersException();
  }
}

template <class Policy> void throwIfNoCoordinates(const DM dm) {
  auto coordinates = Vec{};
  Policy::handleError(DMGetCoordinates(dm, &coordinates));
  if (!coordinates) {
    throw InvalidParametersException();
  }
}
} // namespace detail

template <class Policy>
void writeToViewer(const Mesh<Policy> &mesh, const Viewer<Policy> &viewer) {
  detail::throwIfNoCoordinateDM<Policy>(mesh.data());
  detail::throwIfNoCoordinates<Policy>(mesh.data());
  Policy::handleError(DMView(mesh.data(), viewer.data()));
}

template <class Policy>
void writeToViewer(const distributed<Vector<Policy>> &data,
                   const Viewer<Policy> &viewer) {
  const auto dm = [&]() {
    auto dm = DM{};
    Policy::handleError(VecGetDM(data.unwrap().data(), &dm));
    return dm;
  }();

  if (!dm) {
    throw InvalidParametersException();
  }
  detail::throwIfNoCoordinateDM<Policy>(dm);
  detail::throwIfNoCoordinates<Policy>(dm);

  Policy::handleError(VecView(data.unwrap().data(), viewer.data()));
}

} // namespace cpppetsc
} // namespace ae108
