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

#include "ae108/solve/AffineTransform.h"

namespace ae108 {
namespace solve {

extern template class AffineTransform<cpppetsc::SequentialComputePolicy>;
extern template class AffineTransform<cpppetsc::ParallelComputePolicy>;

template cpppetsc::distributed<
    typename AffineTransform<cpppetsc::SequentialComputePolicy>::vector_type>
apply(const AffineTransform<cpppetsc::SequentialComputePolicy> &,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::SequentialComputePolicy>::vector_type> &);
template cpppetsc::distributed<
    typename AffineTransform<cpppetsc::ParallelComputePolicy>::vector_type>
apply(const AffineTransform<cpppetsc::ParallelComputePolicy> &,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::ParallelComputePolicy>::vector_type> &);

template void
apply(const AffineTransform<cpppetsc::SequentialComputePolicy> &transform,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::SequentialComputePolicy>::vector_type> &,
      cpppetsc::distributed<typename AffineTransform<
          cpppetsc::SequentialComputePolicy>::vector_type> *const);
template void
apply(const AffineTransform<cpppetsc::ParallelComputePolicy> &transform,
      const cpppetsc::distributed<typename AffineTransform<
          cpppetsc::ParallelComputePolicy>::vector_type> &,
      cpppetsc::distributed<typename AffineTransform<
          cpppetsc::ParallelComputePolicy>::vector_type> *const);

} // namespace solve
} // namespace ae108