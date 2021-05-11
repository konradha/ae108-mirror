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

template <> PetscErrorCode callDestructor(VecScatter *const ptr) noexcept {
  return VecScatterDestroy(ptr);
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
} // namespace detail
} // namespace cpppetsc
} // namespace ae108
