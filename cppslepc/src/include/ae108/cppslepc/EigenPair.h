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

#pragma once

#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cppslepc {

template <class Policy> struct EigenPair {
  using vector_type = cpppetsc::distributed<cpppetsc::Vector<Policy>>;
  using value_type = typename vector_type::value_type::value_type;

#ifdef AE108_PETSC_COMPLEX
  value_type value;
  vector_type vector;
#else
  std::complex<value_type> value;
  vector_type vector_real;
  vector_type vector_imag;
#endif
};

} // namespace cppslepc
} // namespace ae108