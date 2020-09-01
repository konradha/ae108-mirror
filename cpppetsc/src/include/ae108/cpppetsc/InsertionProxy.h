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