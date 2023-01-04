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

#include "ae108/cpppetsc/UniqueEntity.h"
#include <petscviewer.h>
#include <petscviewerhdf5.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class Viewer {
public:
  /**
   * @brief Creates a PetscViewer that writes to stdout.
   */
  static Viewer fromStdout();

  /**
   * @brief Used to specify the mode a file is opened with.
   */
  enum class Mode {
    read,
    write,
  };

  /**
   * @brief Creates a PetscViewer that writes to an ASCII file at path.
   */
  static Viewer fromAsciiFilePath(const char *path, const Mode mode);

  /**
   * @brief Creates a PetscViewer that writes to an HDF5 file at path.
   */
  static Viewer fromHdf5FilePath(const char *path, const Mode mode);

  /**
   * @brief Creates a PetscViewer that writes to a Binary file at path.
   */
  static Viewer fromBinaryFilePath(const char *path, const Mode mode);

  /**
   * @brief Creates a Viewer from the provided viewer (takes ownership).
   */
  explicit Viewer(UniqueEntity<PetscViewer> viewer);

  /**
   * @brief Returns the internal PetscViewer.
   */
  PetscViewer data() const;

private:
  UniqueEntity<PetscViewer> _viewer;
};

} // namespace cpppetsc
} // namespace ae108

namespace ae108 {
namespace cpppetsc {

template <class Policy> Viewer<Policy> Viewer<Policy>::fromStdout() {
  auto viewer = PetscViewer{};
  Policy::handleError(
      PetscViewerASCIIGetStdout(Policy::communicator(), &viewer));
  return Viewer(UniqueEntity<PetscViewer>(viewer, [](PetscViewer) {}));
}

template <class Policy>
Viewer<Policy> Viewer<Policy>::fromAsciiFilePath(const char *path,
                                                 const Mode mode) {
  auto viewer = PetscViewer{};

  if (mode == Mode::write) {
    Policy::handleError(
        PetscViewerASCIIOpen(Policy::communicator(), path, &viewer));
  } else {
    Policy::handleError(PetscViewerCreate(Policy::communicator(), &viewer));
    Policy::handleError(PetscViewerSetType(viewer, PETSCVIEWERASCII));
    Policy::handleError(
        PetscViewerFileSetMode(viewer, PetscFileMode::FILE_MODE_READ));
    Policy::handleError(PetscViewerFileSetName(viewer, path));
  }

  return Viewer(makeUniqueEntity<Policy>(viewer));
}

template <class Policy>
Viewer<Policy> Viewer<Policy>::fromHdf5FilePath(const char *path,
                                                const Mode mode) {
  auto viewer = PetscViewer{};
  Policy::handleError(PetscViewerHDF5Open(Policy::communicator(), path,
                                          mode == Mode::write
                                              ? PetscFileMode::FILE_MODE_WRITE
                                              : PetscFileMode::FILE_MODE_READ,
                                          &viewer));
  Policy::handleError(PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_PETSC));
  return Viewer(makeUniqueEntity<Policy>(viewer));
}

template <class Policy>
Viewer<Policy> Viewer<Policy>::fromBinaryFilePath(const char *path,
                                                  const Mode mode) {
  auto viewer = PetscViewer{};
  Policy::handleError(PetscViewerBinaryOpen(Policy::communicator(), path,
                                            mode == Mode::write
                                                ? PetscFileMode::FILE_MODE_WRITE
                                                : PetscFileMode::FILE_MODE_READ,
                                            &viewer));
  return Viewer(makeUniqueEntity<Policy>(viewer));
}

template <class Policy> PetscViewer Viewer<Policy>::data() const {
  return _viewer.get();
}

template <class Policy>
Viewer<Policy>::Viewer(UniqueEntity<PetscViewer> viewer)
    : _viewer(std::move(viewer)) {}

} // namespace cpppetsc
} // namespace ae108