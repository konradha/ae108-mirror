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

#include <functional>
#include <petscsys.h>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief This class is used to provide a "matrix"-like interface (i.e. access
 * elements via operator()) to the Matrix wrapper. PETSc requires that calls to
 * INSERT and ADD are not mixed; therefore specializations of this class only
 * provide one of these methods.
 */
template <InsertMode Mode> struct InsertionProxy {};

/**
 * @brief This class redirects calls to += to the internal functor.
 */
template <> struct InsertionProxy<ADD_VALUES> {
  void operator+=(const PetscScalar value) { _functor(value); }
  const std::function<void(const PetscScalar)> _functor;
};

/**
 * @brief This class redirects calls to = to the internal functor.
 */
template <> struct InsertionProxy<INSERT_VALUES> {
  void operator=(const PetscScalar value) { _functor(value); }
  const std::function<void(const PetscScalar)> _functor;
};
} // namespace cpppetsc
} // namespace ae108