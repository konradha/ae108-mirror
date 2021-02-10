// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/embedding/compute_jacobian.h"
#include "ae108/elements/integrator/IntegrateTrait.h"
#include "ae108/elements/integrator/IntegratorBase.h"
#include "ae108/elements/quadrature/integrate.h"
#include "ae108/elements/tensor/Tensor.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <Eigen/LU>

namespace ae108 {
namespace elements {
namespace integrator {

template <class Shape_, class Quadrature_>
class IsoparametricIntegrator final
    : public IntegratorBase<typename Shape_::size_type,
                            typename Shape_::value_type, Shape_::size(),
                            Shape_::dimension()> {
public:
  using Embedding = embedding::IsoparametricEmbedding<Shape_>;
  using Quadrature = Quadrature_;

  explicit IsoparametricIntegrator(const Embedding &embedding) noexcept {
    auto dxN_iterator = dxN_.begin();
    auto J_iterator = J_.begin();
    for (const auto &point : Quadrature::data.points) {
      const auto gradients = shape::compute_gradients<Shape_>(point);
      const auto jacobian = embedding::compute_jacobian(embedding, point);

      tensor::as_matrix_of_rows(&*(dxN_iterator++)) =
          tensor::as_matrix_of_rows(&gradients) *
          tensor::as_matrix_of_rows(&jacobian).inverse();
      *(J_iterator++) =
          std::abs(tensor::as_matrix_of_rows(&jacobian).determinant());
    }
  }

  using PreTransform =
      tensor::Tensor<typename IsoparametricIntegrator::value_type,
                     Shape_::size(), Embedding::physical_dimension()>;

  const typename Quadrature::template Collection<PreTransform> &pre() const
      noexcept {
    return dxN_;
  };

  using PostTransform = typename IsoparametricIntegrator::value_type;

  const typename Quadrature::template Collection<PostTransform> &post() const
      noexcept {
    return J_;
  };

private:
  typename Quadrature::template Collection<PreTransform> dxN_;
  typename Quadrature::template Collection<PostTransform> J_;
};

template <class Shape_, class Quadrature_>
struct IntegrateTrait<IsoparametricIntegrator<Shape_, Quadrature_>> {
  /**
   * @tparam DiscretizedFunction Equal to an
   * Integrator::DiscretizedFunction<...>.
   */
  template <class Integrator, class R, class F, class DiscretizedFunction>
  R operator()(const Integrator &integrator, F f, const DiscretizedFunction &u,
               R init) const noexcept {
    return quadrature::integrate<Quadrature_>(
        [&f, &u](const typename Quadrature_::size_type id,
                 const typename Quadrature_::Point &,
                 const typename Integrator::PreTransform &pre,
                 const typename Integrator::PostTransform &post) {
          typename Integrator::template Point<
              std::tuple_size<typename DiscretizedFunction::value_type>::value>
              grad_u;
          tensor::as_matrix_of_rows(&grad_u) =
              tensor::as_matrix_of_rows(&u).transpose() *
              tensor::as_matrix_of_rows(&pre);

          return R(post * f(id, grad_u, pre));
        },
        init, integrator.pre(), integrator.post());
  }
};

} // namespace integrator
} // namespace elements
} // namespace ae108