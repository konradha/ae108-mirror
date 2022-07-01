// © 2021, 2022 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/ComputeEnergyTrait.h"
#include "ae108/elements/ComputeForcesTrait.h"
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/tensor/as_vector.h"

namespace ae108 {
namespace elements {

template <class RealType_> struct BarProperties {
  using real_type = RealType_;

  real_type young_modulus;
  real_type cross_section;
};

/**
 * @brief Computes the stiffness matrix for a bar with the given
 * axis and the given properties.
 * @note see Dennis M. Kochmann, "Computational Mechanics I: Introduction to
 * FEA Analysis", https://mm.ethz.ch/education, pp.77-81
 */
template <std::size_t Dimension_, class RealType_>
Eigen::Matrix<RealType_, 2 * Dimension_, 2 * Dimension_, Eigen::RowMajor>
bar_stiffness_matrix(const tensor::Tensor<RealType_, Dimension_> &axis,
                     const BarProperties<RealType_> &properties) noexcept {
  static_assert(Dimension_ > 0);
  /**
   * @brief Computes the 1D stiffness matrix assuming an axis of unit length.
   */
  constexpr auto stiffness_matrix_1D =
      [](const BarProperties<RealType_> &properties) {
        const auto factor = properties.cross_section * properties.young_modulus;
        auto x = Eigen::Matrix<double, 2, 2>();
        x << factor, -1. * factor, -1 * factor, factor;
        return x;
      };

  /**
   * @brief Computes the rotation matrix assuming an axis of unit length.
   */
  constexpr auto rotation =
      [](const tensor::Tensor<RealType_, Dimension_> &axis) {
        auto x = Eigen::Matrix<RealType_, 2, 2 * Dimension_>::Zero().eval();
        x.template block<1, Dimension_>(0, 0) = tensor::as_vector(&axis);
        x.template block<1, Dimension_>(1, Dimension_) =
            tensor::as_vector(&axis);
        return x;
      };

  const auto K = stiffness_matrix_1D(properties);
  const auto R = rotation(axis);

  return std::pow(tensor::as_vector(&axis).norm(), -3.) *
         (R.transpose() * K * R);
}

/**
 * @tparam Dimension_ The dimension of physical space.
 */
template <std::size_t Dimension_, class ValueType_ = double,
          class RealType_ = double>
struct Bar final
    : ElementBase<Bar<Dimension_, ValueType_, RealType_>, std::size_t,
                  ValueType_, RealType_, 2, Dimension_> {
public:
  /**
   * @param matrix A symmetric stiffness matrix.
   */
  explicit Bar(typename Bar::StiffnessMatrix matrix) noexcept
      : stiffness_matrix_(std::move(matrix)) {}

  const typename Bar::StiffnessMatrix &stiffness_matrix() const {
    return stiffness_matrix_;
  }

  /**
   * @brief The dimension of physical space.
   */
  static constexpr typename Bar::size_type dimension() { return Dimension_; }

private:
  typename Bar::StiffnessMatrix stiffness_matrix_;
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeEnergyTrait<Bar<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::Energy
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    const auto v = tensor::as_vector(&u);
    return typename Element::Energy{.5} * v.transpose() *
           element.stiffness_matrix() * v;
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeForcesTrait<Bar<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::Forces
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    typename Element::Forces forces;
    tensor::as_vector(&forces) =
        element.stiffness_matrix() * tensor::as_vector(&u);
    return forces;
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeStiffnessMatrixTrait<Bar<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &element,
             const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return element.stiffness_matrix();
  }
};

} // namespace elements
} // namespace ae108