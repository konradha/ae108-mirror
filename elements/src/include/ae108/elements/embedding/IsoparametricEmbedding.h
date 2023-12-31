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

#include "ae108/elements/embedding/EmbedTrait.h"
#include "ae108/elements/embedding/EmbeddingBase.h"
#include "ae108/elements/embedding/JacobianTrait.h"
#include "ae108/elements/shape/compute_gradients.h"
#include "ae108/elements/shape/compute_values.h"
#include "ae108/elements/tensor/Tensor.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include "ae108/elements/tensor/as_vector.h"
#include <Eigen/Core>
#include <array>
#include <type_traits>

namespace ae108 {
namespace elements {
namespace embedding {

/**
 * @brief An isoparametric embedding based on shape functions: x(xi) = P *
 * N(xi). Here, P is the matrix of points p_j where N_i(p_j) = delta_ij and N is
 * the vector of shape functions.
 */
template <class Shape>
class IsoparametricEmbedding
    : public EmbeddingBase<
          Shape::dimension(), typename Shape::Point, Shape::dimension(),
          tensor::Tensor<typename Shape::value_type, Shape::dimension()>,
          typename Shape::size_type, typename Shape::value_type> {
public:
  template <class T> using Collection = std::array<T, Shape::size()>;

  /**
   * @brief Constructs an embedding from a collection of points in physical
   * space.
   */
  explicit constexpr IsoparametricEmbedding(
      Collection<typename IsoparametricEmbedding::PhysicalPoint>
          points) noexcept
      : points_(std::move(points)) {}

  /**
   * @brief Returns the points in physical space that define the embedding.
   */
  constexpr const Collection<typename IsoparametricEmbedding::PhysicalPoint> &
  points() const noexcept {
    return points_;
  }

private:
  Collection<typename IsoparametricEmbedding::PhysicalPoint> points_;
};

template <class Shape> struct EmbedTrait<IsoparametricEmbedding<Shape>> {
  template <class Embedding>
  typename Embedding::PhysicalPoint
  operator()(const Embedding &embedding,
             const typename Embedding::ReferencePoint &point) const noexcept {
    const auto shape_function_values = shape::compute_values<Shape>(point);
    auto result = typename Embedding::PhysicalPoint();
    const auto &points = embedding.points();

    tensor::as_vector(&result) =
        tensor::as_matrix_of_rows(&points).transpose() *
        tensor::as_vector(&shape_function_values);

    return result;
  }
};

template <class Shape> struct JacobianTrait<IsoparametricEmbedding<Shape>> {
  template <class Embedding>
  using Jacobian = tensor::Tensor<typename Embedding::value_type,
                                  Embedding::physical_dimension(),
                                  Embedding::reference_dimension()>;

  template <class Embedding>
  Jacobian<Embedding>
  operator()(const Embedding &embedding,
             const typename Embedding::ReferencePoint &point) const noexcept {
    const auto shape_function_gradients =
        shape::compute_gradients<Shape>(point);
    const auto &points = embedding.points();

    using ReturnType = tensor::Tensor<typename Embedding::value_type,
                                      Embedding::reference_dimension(),
                                      Embedding::physical_dimension()>;
    auto result = ReturnType();

    tensor::as_matrix_of_rows(&result) =
        tensor::as_matrix_of_rows(&points).transpose() *
        tensor::as_matrix_of_rows(&shape_function_gradients);

    return result;
  }
};

} // namespace embedding
} // namespace elements
} // namespace ae108