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

#include "ae108/elements/ComputeEnergyTrait.h"
#include "ae108/elements/ComputeForcesTrait.h"
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/integrator/integrate.h"
#include "ae108/elements/materialmodels/compute_energy.h"
#include "ae108/elements/materialmodels/compute_stress.h"
#include "ae108/elements/materialmodels/compute_tangent_matrix.h"
#include "ae108/elements/tensor/as_matrix_of_columns.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"

namespace ae108 {
namespace elements {

/**
 * @brief An element that computes energy, forces, and stiffness matrix using a
 * material model and an integrator.
 */
template <class MaterialModel_, class Integrator_>
class CoreElement final
    : public ElementBase<CoreElement<MaterialModel_, Integrator_>,
                         typename Integrator_::size_type,
                         typename Integrator_::value_type, Integrator_::size(),
                         MaterialModel_::degrees_of_freedom(),
                         MaterialModel_::dimension()> {
public:
  using MaterialModel = MaterialModel_;
  using Integrator = Integrator_;

  /**
   * @brief Constructs an element using a material model and an integrator.
   */
  explicit CoreElement(MaterialModel model, Integrator integrator) noexcept
      : model_(std::move(model)), integrator_(std::move(integrator)) {}

  /**
   * @brief The instance of the model used to compute entities at integration
   * points (e.g. energy).
   */
  const MaterialModel &model() const noexcept { return model_; }

  /**
   * @brief The instance of the class used to integrate entities (e.g. energy).
   */
  const Integrator &integrator() const noexcept { return integrator_; }

private:
  MaterialModel model_;
  Integrator integrator_;
};

template <class MaterialModel_, class Integrator_>
struct ComputeEnergyTrait<CoreElement<MaterialModel_, Integrator_>> {
  template <class Element>
  typename Element::Energy
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &time) const noexcept {
    const auto &model = element.model();

    return integrator::integrate<Integrator_>(
        element.integrator(),
        [time,
         &model](const typename Element::size_type id,
                 typename Element::MaterialModel::DisplacementGradient &grad_u,
                 const typename Element::Integrator::PreTransform &) {
          return typename Element::Energy{
              materialmodels::compute_energy(model, id, grad_u, time)};
        },
        u, typename Element::Energy{0.});
  }
};

template <class MaterialModel_, class Integrator_>
struct ComputeForcesTrait<CoreElement<MaterialModel_, Integrator_>> {
  template <class Element>
  typename Element::Forces
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &time) const noexcept {
    auto result = typename Element::Forces();

    const auto &model = element.model();

    auto result_matrix = tensor::as_matrix_of_columns(&result);
    result_matrix = integrator::integrate<Integrator_>(
        element.integrator(),
        [time,
         &model](const typename Element::size_type id,
                 typename Element::MaterialModel::DisplacementGradient &grad_u,
                 const typename Element::Integrator::PreTransform &dxN) {
          const auto stress =
              materialmodels::compute_stress(model, id, grad_u, time);
          return (tensor::as_matrix_of_rows(&stress) *
                  tensor::as_matrix_of_rows(&dxN).transpose())
              .eval();
        },
        u, result_matrix.Zero().eval());

    return result;
  }
};

template <class MaterialModel_, class Integrator_>
struct ComputeStiffnessMatrixTrait<CoreElement<MaterialModel_, Integrator_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &time) const noexcept {
    using size_type = typename Element::size_type;

    const auto &model = element.model();

    return integrator::integrate<Integrator_>(
        element.integrator(),
        [time,
         &model](const typename Element::size_type id,
                 typename Element::MaterialModel::DisplacementGradient &grad_u,
                 const typename Element::Integrator::PreTransform &dxN) {
          constexpr auto degrees_of_freedom = Element::degrees_of_freedom();

          auto result = Element::StiffnessMatrix::Zero().eval();
          const auto tangent =
              materialmodels::compute_tangent_matrix(model, id, grad_u, time);

          for (auto i = size_type{0}; i < degrees_of_freedom; ++i) {
            for (auto k = size_type{0}; k < degrees_of_freedom; ++k) {
              auto result_slice = ResultSlice<Element>(&result(i, k));
              auto tangent_slice = TangentSlice<Element>(
                  degrees_of_freedom > 0 ? tangent[i].front()[k].data()
                                         : nullptr);

              result_slice = tensor::as_matrix_of_rows(&dxN) * tangent_slice *
                             tensor::as_matrix_of_rows(&dxN).transpose();
            }
          }
          return result;
        },
        u, Element::StiffnessMatrix::Zero().eval());
  }

private:
  template <class Element>
  using ResultSlice =
      Eigen::Map<Eigen::Matrix<typename Element::value_type, Element::size(),
                               Element::size(), Eigen::RowMajor>,
                 0,
                 Eigen::Stride<Eigen::Index{Element::size() *
                                            Element::degrees_of_freedom() *
                                            Element::degrees_of_freedom()},
                               Eigen::Index{Element::degrees_of_freedom()}>>;

  template <class Element>
  using TangentSlice = typename std::enable_if<
      Element::degrees_of_freedom() == 0 ||
          sizeof(typename Element::MaterialModel::TangentMatrix) ==
              Element::degrees_of_freedom() * Element::dimension() *
                  Element::degrees_of_freedom() * Element::dimension() *
                  sizeof(typename Element::MaterialModel::value_type),
      Eigen::Map<
          const Eigen::Matrix<typename Element::MaterialModel::value_type,
                              Element::dimension(), Element::dimension(),
                              Eigen::RowMajor>,
          0,
          Eigen::Stride<Eigen::Index{Element::dimension() *
                                     Element::degrees_of_freedom()},
                        Eigen::Index{1}>>>::type;
};

} // namespace elements
} // namespace ae108