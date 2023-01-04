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

#include "ae108/cpppetsc/UniqueEntity.h"
#include <petscdm.h>
#include <petscis.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsection.h>
#include <petscsf.h>
#include <petscsnes.h>
#include <petsctao.h>
#include <petscvec.h>
#include <petscviewer.h>

namespace ae108 {
namespace cpppetsc {
namespace detail {

template <> PetscErrorCode callDestructor(IS *const ptr) noexcept {
  return ISDestroy(ptr);
}

template <> PetscErrorCode callDestructor(Vec *const ptr) noexcept {
  return VecDestroy(ptr);
}

template <> PetscErrorCode callDestructor(PetscSF *const ptr) noexcept {
  return PetscSFDestroy(ptr);
}

template <> PetscErrorCode callDestructor(PetscSection *const ptr) noexcept {
  return PetscSectionDestroy(ptr);
}

template <> PetscErrorCode callDestructor(DM *const ptr) noexcept {
  return DMDestroy(ptr);
}

template <> PetscErrorCode callDestructor(Mat *const ptr) noexcept {
  return MatDestroy(ptr);
}

template <> PetscErrorCode callDestructor(SNES *const ptr) noexcept {
  return SNESDestroy(ptr);
}

template <> PetscErrorCode callDestructor(KSP *const ptr) noexcept {
  return KSPDestroy(ptr);
}

template <> PetscErrorCode callDestructor(Tao *const ptr) noexcept {
  return TaoDestroy(ptr);
}

template <> PetscErrorCode callDestructor(PetscViewer *const ptr) noexcept {
  return PetscViewerDestroy(ptr);
}

template <>
PetscErrorCode callDestructor(ISLocalToGlobalMapping *const ptr) noexcept {
  return ISLocalToGlobalMappingDestroy(ptr);
}
} // namespace detail
} // namespace cpppetsc
} // namespace ae108
