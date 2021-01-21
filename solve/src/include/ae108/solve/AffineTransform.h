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