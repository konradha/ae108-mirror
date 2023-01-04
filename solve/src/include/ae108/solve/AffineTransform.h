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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/clone.h"

namespace ae108 {
namespace solve {

template <class Policy> struct AffineTransform {
  using matrix_type = cpppetsc::Matrix<Policy>;
  using vector_type = cpppetsc::Vector<Policy>;

  matrix_type matrix;
  cpppetsc::distributed<vector_type> vector;
};

/**
 * @brief Applies the affine transform to an input vector `x` and returns the
 * result.
 */
template <class Policy>
cpppetsc::distributed<typename AffineTransform<Policy>::vector_type>
apply(const AffineTransform<Policy> &transform,
      const cpppetsc::distributed<typename AffineTransform<Policy>::vector_type>
          &x) {
  auto result = clone(transform.vector);
  result.unwrap().addAx(transform.matrix, x);
  return result;
}

/**
 * @brief Applies the affine transform to an input vector `x` and writes the
 * result to `out`.
 */
template <class Policy>
void apply(
    const AffineTransform<Policy> &transform,
    const cpppetsc::distributed<typename AffineTransform<Policy>::vector_type>
        &x,
    cpppetsc::distributed<typename AffineTransform<Policy>::vector_type>
        *const out) {
  copy(transform.vector, out);
  out->unwrap().addAx(transform.matrix, x);
}

extern template class AffineTransform<cpppetsc::SequentialComputePolicy>;
extern template class AffineTransform<cpppetsc::ParallelComputePolicy>;

extern template cpppetsc::distributed<
    typename AffineTransform<cpppetsc::SequentialComputePolicy>::vector_type>
apply(const AffineTransform<cpppetsc::SequentialComputePolicy> &,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::SequentialComputePolicy>::vector_type> &);
extern template cpppetsc::distributed<
    typename AffineTransform<cpppetsc::ParallelComputePolicy>::vector_type>
apply(const AffineTransform<cpppetsc::ParallelComputePolicy> &,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::ParallelComputePolicy>::vector_type> &);

extern template void
apply(const AffineTransform<cpppetsc::SequentialComputePolicy> &transform,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::SequentialComputePolicy>::vector_type> &,
      cpppetsc::distributed<typename AffineTransform<
          cpppetsc::SequentialComputePolicy>::vector_type> *const);
extern template void
apply(const AffineTransform<cpppetsc::ParallelComputePolicy> &transform,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::ParallelComputePolicy>::vector_type> &,
      cpppetsc::distributed<typename AffineTransform<
          cpppetsc::ParallelComputePolicy>::vector_type> *const);

} // namespace solve
} // namespace ae108