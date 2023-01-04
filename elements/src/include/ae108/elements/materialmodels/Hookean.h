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

#include "ae108/elements/materialmodels/ComputeEnergyTrait.h"
#include "ae108/elements/materialmodels/ComputeStrainTrait.h"
#include "ae108/elements/materialmodels/ComputeStressTrait.h"
#include "ae108/elements/materialmodels/ComputeTangentMatrixTrait.h"
#include "ae108/elements/materialmodels/MaterialModelBase.h"
#include "ae108/elements/materialmodels/compute_strain.h"
#include "ae108/elements/materialmodels/compute_stress.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"

namespace ae108 {
namespace elements {
namespace materialmodels {

/**
 * @brief A material model for linear elastic isotropic materials.
 * @remark See https://en.wikipedia.org/wiki/Hooke%27s_law
 */
template <std::size_t Dimension_, class ValueType_ = double,
          class RealType_ = double>
class Hookean final
    : public MaterialModelBase<std::size_t, ValueType_, RealType_, Dimension_> {
public:
  /**
   * @brief Constructs the material model using Young's modulus E and the
   * Poisson ratio nu.
   */
  explicit constexpr Hookean(
      const typename Hookean::real_type &youngs_modulus,
      const typename Hookean::real_type &poisson_ratio) noexcept
      : lambda_(youngs_modulus * poisson_ratio /
                (typename Hookean::real_type{1.} -
                 typename Hookean::real_type{2.} * poisson_ratio) /
                (typename Hookean::real_type{1.} + poisson_ratio)),
        mu_(youngs_modulus / typename Hookean::real_type{2.} /
            (typename Hookean::real_type{1.} + poisson_ratio)) {}

  /**
   * @brief Returns the first Lame constant.
   */
  constexpr typename Hookean::real_type lambda() const noexcept {
    return lambda_;
  }

  /**
   * @brief Returns the second Lame constant.
   */
  constexpr typename Hookean::real_type mu() const noexcept { return mu_; }

private:
  typename Hookean::real_type lambda_;
  typename Hookean::real_type mu_;
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeEnergyTrait<Hookean<Dimension_, ValueType_, RealType_>> {
  template <class MaterialModel>
  typename MaterialModel::Energy
  operator()(const MaterialModel &model,
             const typename MaterialModel::size_type id,
             const typename MaterialModel::DisplacementGradient &gradient,
             const typename MaterialModel::Time time) noexcept {
    const auto strain = compute_strain(model, id, gradient, time);
    const auto stress = compute_stress(model, id, gradient, time);
    return std::real(typename MaterialModel::real_type{0.5} *
                     tensor::as_matrix_of_rows(&strain)
                         .cwiseProduct(tensor::as_matrix_of_rows(&stress))
                         .sum());
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeStrainTrait<Hookean<Dimension_, ValueType_, RealType_>> {
  template <class MaterialModel>
  typename MaterialModel::Strain
  operator()(const MaterialModel &, const typename MaterialModel::size_type,
             const typename MaterialModel::DisplacementGradient &gradient,
             const typename MaterialModel::Time) const noexcept {
    const auto gradient_matrix = tensor::as_matrix_of_rows(&gradient);

    typename MaterialModel::Strain result;
    tensor::as_matrix_of_rows(&result) =
        typename MaterialModel::real_type{.5} *
        (gradient_matrix + gradient_matrix.adjoint());
    return result;
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeStressTrait<Hookean<Dimension_, ValueType_, RealType_>> {
  template <class MaterialModel>
  typename MaterialModel::Stress
  operator()(const MaterialModel &model,
             const typename MaterialModel::size_type id,
             const typename MaterialModel::DisplacementGradient &gradient,
             const typename MaterialModel::Time time) const noexcept {
    const auto strain = compute_strain(model, id, gradient, time);
    const auto strain_matrix = tensor::as_matrix_of_rows(&strain);

    typename MaterialModel::Stress result;
    tensor::as_matrix_of_rows(&result) =
        model.lambda() * strain_matrix.trace() * strain_matrix.Identity() +
        typename MaterialModel::real_type{2.} * model.mu() * strain_matrix;
    return result;
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeTangentMatrixTrait<Hookean<Dimension_, ValueType_, RealType_>> {
  template <class MaterialModel>
  typename MaterialModel::TangentMatrix
  operator()(const MaterialModel &model,
             const typename MaterialModel::size_type,
             const typename MaterialModel::DisplacementGradient &,
             const typename MaterialModel::Time) const noexcept {
    using size_type = typename MaterialModel::size_type;

    typename MaterialModel::TangentMatrix result;

    for (auto i = size_type{0}; i < MaterialModel::dimension(); ++i) {
      for (auto j = size_type{0}; j < MaterialModel::dimension(); ++j) {
        auto matrix = tensor::as_matrix_of_rows(&result[i][j]);
        matrix = (i == j ? (model.lambda() * matrix.Identity()).eval()
                         : matrix.Zero().eval());
        matrix(i, j) += model.mu();
        matrix(j, i) += model.mu();
      }
    }

    return result;
  }
};

} // namespace materialmodels
} // namespace elements
} // namespace ae108
