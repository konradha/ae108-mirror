// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/embedding/compute_jacobian.h"
#include "ae108/elements/integrator/ComputeVolumeTrait.h"
#include "ae108/elements/integrator/IntegrateShapeTrait.h"
#include "ae108/elements/integrator/IntegrateTrait.h"
#include "ae108/elements/integrator/IntegratorBase.h"
#include "ae108/elements/quadrature/integrate.h"
#include "ae108/elements/tensor/Tensor.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <Eigen/LU>
#include <algorithm>

namespace ae108 {
namespace elements {
namespace integrator {

template <class Shape_, class Quadrature_, class ValueType_ = double,
          class RealType_ = double>
class IsoparametricIntegrator final
    : public IntegratorBase<typename Shape_::size_type, ValueType_, RealType_,
                            Shape_::size(), Shape_::dimension()> {
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

  const typename Quadrature::template Collection<PreTransform> &
  pre() const noexcept {
    return dxN_;
  };

  using PostTransform = typename IsoparametricIntegrator::real_type;

  const typename Quadrature::template Collection<PostTransform> &
  post() const noexcept {
    return J_;
  };

private:
  typename Quadrature::template Collection<PreTransform> dxN_;
  typename Quadrature::template Collection<PostTransform> J_;
};

template <class Shape_, class Quadrature_, class ValueType_, class RealType_>
struct IntegrateTrait<
    IsoparametricIntegrator<Shape_, Quadrature_, ValueType_, RealType_>> {
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

template <class Shape_, class Quadrature_, class ValueType_, class RealType_>
struct IntegrateShapeTrait<
    IsoparametricIntegrator<Shape_, Quadrature_, ValueType_, RealType_>> {
  template <class Integrator, class R, class F>
  R operator()(const Integrator &integrator, F f, R init) const noexcept {
    using Value =
        typename Shape_::template Collection<typename Shape_::value_type>;
    auto values = typename Quadrature_::template Collection<Value>();

    std::transform(Quadrature_::data.points.begin(),
                   Quadrature_::data.points.end(), values.begin(),
                   &shape::compute_values<Shape_>);

    return quadrature::integrate<Quadrature_>(
        [&f](const typename Quadrature_::size_type id,
             const typename Quadrature_::Point &,
             const typename Integrator::PostTransform &post,
             const Value &value) { return R(post * f(id, value)); },
        init, integrator.post(), values);
  }
};

template <class Shape_, class Quadrature_, class ValueType_, class RealType_>
struct ComputeVolumeTrait<
    IsoparametricIntegrator<Shape_, Quadrature_, ValueType_, RealType_>> {
  template <class Integrator>
  auto operator()(const Integrator &integrator) const noexcept {
    return quadrature::integrate<Quadrature_>(
        [](auto &&, auto &&, const auto &post) { return post; },
        typename Integrator::real_type{0.}, integrator.post());
  }
};

} // namespace integrator
} // namespace elements
} // namespace ae108
