// © 2023 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/ComputeEnergyTrait.h"
#include "ae108/elements/ComputeForcesTrait.h"
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/tensor/as_vector.h"

#include "ae108/elements/materialmodels/ComputeStressTrait.h"
#include "ae108/elements/materialmodels/ComputeTangentMatrixTrait.h"

// #include "ae108/elements/TwoNodeCorotationalBeamJHElement.h"
// #include "ae108/elements/tensor/as_matrix_of_rows.h"


namespace ae108 {
namespace elements {

template <class RealType_, std::size_t Dimension_>
struct TwoNodeCorotationalBeamJHProperties;

template <class RealType_>
struct TwoNodeCorotationalBeamJHProperties<RealType_, 2> {
  using real_type = RealType_;

  real_type modulus;
  real_type areaMoment;
  real_type area;
  real_type density;
  tensor::Tensor<double, 4> local_angles;
  // real_type torsion_x; leave for later
  // real_type torsion_y;
};


template <std::size_t Dimension_, class Element>
tensor::Tensor<double, Dimension_ *(Dimension_ + 1),
               Dimension_ *(Dimension_ + 1)>
stiffness_matrix(const TwoNodeCorotationalBeamJHProperties<double, Dimension_> &properties,
                 const typename Element::NodalDisplacements &displacements) noexcept;

template <std::size_t Dimension_, class Element>
tensor::Tensor<double, 6, 6>
stiffness_matrix(const TwoNodeCorotationalBeamJHProperties<double, 2> &properties,
                    const typename Element::NodalDisplacements &displacements,
                    const double length) noexcept {
  auto tangent = Element::ComputeTangentMatrixTrait<MaterialModel_>(model, id, gradient); // TODO
  double cosbeta, sinbeta;
  cosbeta = properties.local_angles[0]; sinbeta = properties.local_angles[1];
  auto stress = Element::ComputeStressTrait<Dimension_>();


  auto L  = length;
  // assume all matrices local here; no left/right stress; only entire StressTrait

  auto rg = properties.areaMoment / properties.area;
  auto n  = properties.area * stress[0];
  auto k1 = tangent[0, 0] * properties.area;

  auto c1 = tangent[0, 0] * properties.area;
  auto c2 = k1 * 2. * rg; auto c3 = 2. * c2;
  auto _ = 0;
  tensor::Tensor<double, 3, 3> C = 
  {{
      {{  c1,  c3, _}},
      {{  c3,  c2,  _}},
      {{  _,   _,  _}}
  }};
  // Dimenion_ == DOF? number of nodes?
  tensor::Tensor<double, 4, 1> r = 
  {{
    {{sinbeta, -cosbeta, -sinbeta, cosbeta}}
  }};

  tensor::Tensor<double, 4, 1> z =
  {{
    {{-cosbeta, -sinbeta, cosbeta, sinbeta}}
  }};

  tensor::Tensor<double, Dimension_*(Dimension_+1), Dimension_*(Dimension_+1)> k =
  C + 1./L * z.transpose() * z + ((m1 + m2) / L / L) * (r * z.tranpose()) + z * r.transpose(); 

  return k;
}


template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
twonode_corotational_beamjh_stiffness_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TwoNodeCorotationalBeamJHProperties<double, Dimension_> &properties) noexcept {
      Eigen::Matrix<double, Dimension_ *(Dimension_ + 1), Dimension_ *(Dimension_ + 1), Eigen::RowMajor> values;
      return values;
    }


template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
twonode_corotational_beamjh_lumped_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TwoNodeCorotationalBeamJHProperties<double, Dimension_> &properties,
    const double density) noexcept;


template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
twonode_corotational_beamjh_consistent_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TwoNodeCorotationalBeamJHProperties<double, Dimension_> &properties,
    const double density) noexcept;


template <std::size_t Dimension_, class ValueType_ = double,
          class RealType_ = double>
struct TwoNodeCorotationalBeamJHElement final
    : ElementBase<TwoNodeCorotationalBeamJHElement<Dimension_, ValueType_, RealType_>,
                  std::size_t, ValueType_, RealType_, 2, Dimension_,
                  (Dimension_ * (Dimension_ + 1)) / 2> {
public:
  explicit TwoNodeCorotationalBeamJHElement(
      typename TwoNodeCorotationalBeamJHElement::StiffnessMatrix matrix) noexcept
      : stiffness_matrix_(std::move(matrix)) {}

  const typename TwoNodeCorotationalBeamJHElement::StiffnessMatrix &
  stiffness_matrix() const {
    return stiffness_matrix_;
  }

private:
  typename TwoNodeCorotationalBeamJHElement::StiffnessMatrix stiffness_matrix_;
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeEnergyTrait<
    TwoNodeCorotationalBeamJHElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::Energy
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    const auto v = tensor::as_vector(&u);
    return typename Element::Energy{.5} * v.transpose() *
           element.stiffness_matrix() * v; // placeholder
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeForcesTrait<
    TwoNodeCorotationalBeamJHElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::Forces
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    typename Element::Forces forces;
    tensor::as_vector(&forces) =
        element.stiffness_matrix() * tensor::as_vector(&u); // placeholder
    return forces;
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeStiffnessMatrixTrait<
    TwoNodeCorotationalBeamJHElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &element,
             const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return element.stiffness_matrix();
  }
};



template <class RealType_>
struct TwoNodeCorotationalBeamJHProperties<RealType_, 3> {
  using real_type = RealType_;
};



}
}