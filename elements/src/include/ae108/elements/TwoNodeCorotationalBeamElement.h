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

// #include "ae108/elements/TwoNodeCorotationalBeamElement.h"
// #include "ae108/elements/tensor/as_matrix_of_rows.h"


namespace ae108 {
namespace elements {

template <class RealType_, std::size_t Dimension_>
struct TwoNodeCorotationalBeamProperties;

template <class RealType_>
struct TwoNodeCorotationalBeamProperties<RealType_, 2> {
  using real_type = RealType_;

  real_type modulus; // TODO find out `nature` of modulus
  real_type area_moment;
  real_type area;
  real_type density;


  // double _modulus;
  // double _areaMoment;
  // double _area;
  // double _density;

  // ↓↓↓ COMPARISON ↓↓↓

  // CURRENT Timoshenko element
  // real_type young_modulus;
  // real_type shear_modulus;
  // real_type shear_correction_factor_y;
  // real_type area;
  // real_type area_moment_z;

  // LEGACY Timoshenko element
  // double _shearcorrection;
  // double _area;
  // double _areaMomenty;
  // double _areaMomentz;
  // double _JMoment;
  // Matrix<double,3,1> _desiredLocalZOrientation;
  // double _maxExtensiony;
  // double _maxExtensionz;
  // double _yieldStress;
};

template <class RealType_, class Element>
tensor::Tensor<RealType_, 2> computeLocAxialValues(
  const typename Element::NodalDisplacements &s, // TODO: needs to be changed
  const typename Element::NodalDisplacements &displacements) noexcept
{
  tensor::Tensor<RealType_, 2> x0, x1;
  x0 = s(0) + displacements(0); // TODO: check if this makes sense in matmul backend notation
  x1 = s(1) + displacements(1); // TODO: check if this makes sense in matmul backend notation
  auto d = tensor::as_vector(x1 - x0);
  RealType_ dist = d.norm();
  
  RealType_ axialDisplacement = (dist*dist - x0*x0) / (dist + x0); // what are we actually calculating here?
  return {axialDisplacement, dist};


  /*
  // COMPUTE L
  Point defX0;
  defX0(0, 0) = _x0(0) + displacements[0](0);
  defX0(1, 0) = _x0(1) + displacements[0](1);
  Point defX1;
  defX1(0, 0) = _x1(0) + displacements[1](0);
  defX1(1, 0) = _x1(1) + displacements[1](1);
  Point defDistanceVector = defX1 - defX0;
  double defDistance = defDistanceVector.norm();
  // COMPUTE UL
  double locAxialDisplacement =
      (pow(defDistance, 2) - pow(_distance, 2)) / (defDistance + _distance);
  Vector2d locAxialValues;
  locAxialValues(0) = locAxialDisplacement;
  locAxialValues(1) = defDistance;
  return locAxialValues;
  */
}


template <std::size_t Dimension_, class Element>
tensor::Tensor<double, Dimension_ *(Dimension_ + 1),
               Dimension_ *(Dimension_ + 1)>
stiffness_matrix(const TwoNodeCorotationalBeamProperties<double, Dimension_> &properties,
                 const typename Element::NodalDisplacements &displacements) noexcept;

template <std::size_t Dimension_, class Element>
tensor::Tensor<double, 12, 12>
stiffness_matrix(const TwoNodeCorotationalBeamProperties<double, 3> &properties,
                    const typename Element::NodalDisplacements &displacements,
                    const double length) noexcept {
  // const double L = length;
  // computeLocAxialValues(displacements);
  // const auto M   = properties.modulus; 
  // const auto I   = properties.area_moment;
  // const auto A   = properties.area;

  // compute_strain()

  // // Calculate local angles, bending moments, axial strain
  // auto localAngles = calculateLocalAngles(displacements);
  // auto moments     = calculateLocalMoments(localAngles, M, I);
  // auto axialValues = computeLocAxialValues(displacements);
  // const auto axialStrain = axialValues[0] / L;
  
  // // assembly final stiffness matrix  
  // tensor::Tensor<double, Dimension_*(Dimension_+1), Dimension_*(Dimension_+1)> k;
  
  // // Axial stiffness portion
  // k(0,0) = M*A / L; 
  // k(0,6) = -k(0,0);

  // // Bending stiffness portion
  // // Use properties, localAngles, and moments

  // return k;
  tensor::Tensor<double, Dimension_*(Dimension_+1), Dimension_*(Dimension_+1)> k;
  return k;
}



// template <std::size_t Dimension_, class ValueType_, class RealType_>
// struct ComputeEnergyTrait<
//     TwoNodeCorotationalBeamElement<Dimension_, ValueType_, RealType_>> {
//   template <class Element>
//   typename Element::Energy
//   operator()(const Element &element,
//              const typename Element::NodalDisplacements &u,
//              const typename Element::Time &) const noexcept {
//     const auto v = tensor::as_vector(&u);
//     return typename Element::Energy{.5} * v.transpose() *
//            element.stiffness_matrix() * v;
//   }
// };


template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
twonode_corotational_beam_stiffness_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TwoNodeCorotationalBeamProperties<double, Dimension_> &properties) noexcept;


template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
twonode_corotational_beam_lumped_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TwoNodeCorotationalBeamProperties<double, Dimension_> &properties,
    const double density) noexcept;


template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
twonode_corotational_beam_consistent_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TwoNodeCorotationalBeamProperties<double, Dimension_> &properties,
    const double density) noexcept;


template <std::size_t Dimension_, class ValueType_ = double,
          class RealType_ = double>
struct TwoNodeCorotationalBeamElement final
    : ElementBase<TwoNodeCorotationalBeamElement<Dimension_, ValueType_, RealType_>,
                  std::size_t, ValueType_, RealType_, 2, Dimension_,
                  (Dimension_ * (Dimension_ + 1)) / 2> {
public:
  explicit TwoNodeCorotationalBeamElement(
      typename TwoNodeCorotationalBeamElement::StiffnessMatrix matrix) noexcept
      : stiffness_matrix_(std::move(matrix)) {}

  const typename TwoNodeCorotationalBeamElement::StiffnessMatrix &
  stiffness_matrix() const {
    return stiffness_matrix_;
  }

private:
  typename TwoNodeCorotationalBeamElement::StiffnessMatrix stiffness_matrix_;
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeEnergyTrait<
    TwoNodeCorotationalBeamElement<Dimension_, ValueType_, RealType_>> {
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
    TwoNodeCorotationalBeamElement<Dimension_, ValueType_, RealType_>> {
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
    TwoNodeCorotationalBeamElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &element,
             const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return element.stiffness_matrix();
  }
};



template <class RealType_>
struct TwoNodeCorotationalBeamProperties<RealType_, 3> {
  using real_type = RealType_;
};



}
}