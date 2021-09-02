// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include <complex>
#include <petscsys.h>
#include <type_traits>

int main() {
#ifdef AE108_PETSC_COMPLEX
  static_assert(std::is_same<PetscScalar, std::complex<double>>::value, "");
#endif

#ifdef AE108_PETSC_REAL
  static_assert(std::is_same<PetscScalar, double>::value, "");
#endif
}