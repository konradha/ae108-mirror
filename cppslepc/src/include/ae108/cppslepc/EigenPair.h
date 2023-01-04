// © 2021 ETH Zurich, Mechanics and Materials Lab
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

extern template class EigenPair<cpppetsc::SequentialComputePolicy>;
extern template class EigenPair<cpppetsc::ParallelComputePolicy>;
} // namespace cppslepc
} // namespace ae108