// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#include "ae108/cpppetsc/IndexBasedAccessIterator.h"
#include "ae108/cpppetsc/InsertionProxy.h"
#include "ae108/cpppetsc/IteratorRange.h"
#include "ae108/cpppetsc/Matrix_fwd.h"
#include "ae108/cpppetsc/Mesh_fwd.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector_fwd.h"
#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <memory>
#include <numeric>
#include <petscis.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petscvec.h>
#include <utility>
#include <vector>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class Vector {
public:
  using size_type = PetscInt;
  using value_type = PetscScalar;
  using matrix_type = Matrix<Policy>;

  using const_iterator = IndexBasedAccessIterator<Vector>;

  /**
   * @brief Allocates a vector of length global_size.
   *
   * @remark Does not initialize the vector.
   */
  explicit Vector(const size_type global_size);

  /**
   * @brief Allocates a vector of length global_size and initializes the vector
   * to value.
   */
  explicit Vector(const size_type global_size, const value_type value);

  /**
   * @brief Constructs a vector with the provided entries.
   */
  static Vector fromList(const std::initializer_list<value_type> list);

  /**
   * @brief Creates a global vector with a structure defined by the given mesh.
   * The vector is zero-initialized.
   */
  static distributed<Vector> fromGlobalMesh(const Mesh<Policy> &mesh);

  /**
   * @brief Creates a local vector with a structure defined by the given mesh.
   * The vector is zero-initialized.
   */
  static local<Vector> fromLocalMesh(const Mesh<Policy> &mesh);

  /**
   * @brief Creates a new full sequential copy of the vector.
   */
  static global<Vector> fromDistributed(const distributed<Vector> &vector);

  /**
   * @brief Creates a new full sequential copy of the vector.
   *
   * @param vertexDof Specifies the number of degrees of freedom that
   * every vertex has.
   *
   * @param elementDof Specifies the maximum number of degrees of freedom
   * that every element has.
   *
   * @return The returned vector has size vertexDof * totalNumberOfVertices +
   * elementDof * totalNumberOfElements.
   *
   * @pre The vector is in the default "Mesh"-based ordering.
   * @post The vector is in the canonical ordering:
   * [dof 0 for vertex 0, dof 1 for vertex 0, ..., dof 0 for vertex 1, ...,
   * dof 0 for element 0, dof 1 for element 0, ..., dof 0 for element 1, ...]
   */
  static global<Vector> fromDistributedInCanonicalOrder(
      const distributed<Vector> &vector, const Mesh<Policy> &mesh,
      const size_type vertexDof, const size_type elementDof = 0);

  /**
   * @brief Returns a zero initialized vector with the same layout as the
   * parameter.
   */
  static Vector fromLayoutOf(const Vector &vector);

  /**
   * @brief Returns a deep copy of the vector.
   */
  static Vector clone(const Vector &vector);

  /**
   * @brief Wraps the given UniqueEntity in a Vector.
   */
  explicit Vector(UniqueEntity<Vec> vec);

  /**
   * @brief Provides element() and elements() access functions to set/add
   * single/multiple elements. These calls can be chained (vector is a Vector
   * object)
   *
   * \code{.cpp}vector.replace().element(0, 7.).element(1, 3.);\endcode
   *
   * It is also possible to set/add entries by accessing the index first:
   *
   * \code{.cpp}instance.replace()(0) = 7.;\endcode
   *
   * This way, the Inserter class can be used as a "vector" for assembly.
   */
  template <InsertMode Mode> class Inserter;

  /**
   * @brief Replace elements with different elements. See the documentation of
   * the Inserter class for more information on how to use the return value.
   *
   * @remark Automatically broadcasts changes. Do not use add and replace
   * Inserters at the same time.
   */
  Inserter<INSERT_VALUES> replace();

  /**
   * @brief Add elements to current elements. See the documentation of the
   * Inserter class for more information on how to use the return value.
   *
   * @remark Automatically broadcasts changes. Do not use add and replace
   * Inserters at the same time.
   */
  Inserter<ADD_VALUES> add();

  /**
   * @brief Computes the 2-norm of the vector.
   */
  value_type norm() const;

  /**
   * @brief Print the vector to world stdout.
   */
  void print() const;

  /**
   * @brief Access the element at index.
   */
  value_type operator()(const size_type index) const;

  /**
   * @brief Access the element at index.
   */
  value_type operator[](const size_type index) const;

  /**
   * @brief Scale the vector by factor.
   */
  void scale(const value_type factor);

  /**
   * @brief Sets the vector to alpha * (current value) + beta * x.
   */
  void timesAlphaPlusBetaX(const value_type alpha, const value_type beta,
                           const Vector &x);

  /**
   * @brief Sets the vector to alpha * (current value) + beta * x + gamma * y.
   */
  void timesAlphaPlusBetaXPlusGammaY(const value_type alpha,
                                     const value_type beta, const Vector &x,
                                     const value_type gamma, const Vector &y);

  /**
   * @brief Adds A * x to the vector.
   */
  void addAx(const matrix_type &A, const distributed<Vector> &x);

  /**
   * @brief Returns the global size of the vector.
   */
  size_type size() const;

  /**
   * @brief Returns the local row indices in the format:
   *
   * (first row index, last row index + 1).
   */
  std::pair<size_type, size_type> localRowRange() const;

  /**
   * @brief The begin of the local row range.
   */
  const_iterator localBegin() const;

  /**
   * @brief The end of the local row range.
   */
  const_iterator localEnd() const;

  /**
   * @brief Returns a loopable range of the local values.
   */
  IteratorRange<const_iterator> localRange() const;

  /**
   * @brief Set all elements of the vector to value.
   */
  void fill(const value_type value);

  /**
   * @brief Replaces all entries by zero.
   */
  void setZero();

  /**
   * @brief Returns the internal vector.
   *
   * @remark For internal use only.
   */
  Vec data() const;

private:
  static Vec createVec(const size_type global_size);

  UniqueEntity<Vec> _vec;
};

extern template class Vector<SequentialComputePolicy>;
extern template class Vector<ParallelComputePolicy>;
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/Mesh.h"
#include <petscdm.h>
#include <petscmat.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy>
Vec Vector<Policy>::createVec(const size_type global_size) {
  Vec vec;
  Policy::handleError(
      VecCreateMPI(Policy::communicator(), PETSC_DECIDE, global_size, &vec));
  return vec;
}

template <class Policy>
Vector<Policy>::Vector(const size_type global_size)
    : _vec(makeUniqueEntity<Policy>(createVec(global_size))) {
  Policy::handleError(VecSetFromOptions(_vec.get()));
}

template <class Policy>
Vector<Policy>::Vector(const size_type global_size, const value_type value)
    : Vector(global_size) {
  fill(value);
}

template <class Policy>
Vector<Policy>::Vector(UniqueEntity<Vec> vec) : _vec(std::move(vec)) {}

template <class Policy>
distributed<Vector<Policy>>
Vector<Policy>::fromGlobalMesh(const Mesh<Policy> &mesh) {
  auto vec = Vec();
  Policy::handleError(DMCreateGlobalVector(mesh.data(), &vec));
  auto vector = Vector(makeUniqueEntity<Policy>(vec));

  vector.setZero();

  return tag<DistributedTag>(std::move(vector));
}

template <class Policy>
local<Vector<Policy>> Vector<Policy>::fromLocalMesh(const Mesh<Policy> &mesh) {
  auto vec = Vec();
  Policy::handleError(DMCreateLocalVector(mesh.data(), &vec));
  auto vector = Vector(makeUniqueEntity<Policy>(vec));

  vector.setZero();

  return tag<LocalTag>(std::move(vector));
}

template <class Policy>
global<Vector<Policy>>
Vector<Policy>::fromDistributed(const distributed<Vector> &vector) {
  auto vec = Vec();
  auto vecScatter = VecScatter();
  Policy::handleError(
      VecScatterCreateToAll(vector.unwrap()._vec.get(), &vecScatter, &vec));
  auto fullVector = Vector(makeUniqueEntity<Policy>(vec));
  auto scatterPtr = makeUniqueEntity<Policy>(vecScatter);
  Policy::handleError(
      VecScatterBegin(scatterPtr.get(), vector.unwrap()._vec.get(),
                      fullVector._vec.get(), INSERT_VALUES, SCATTER_FORWARD));
  Policy::handleError(
      VecScatterEnd(scatterPtr.get(), vector.unwrap()._vec.get(),
                    fullVector._vec.get(), INSERT_VALUES, SCATTER_FORWARD));
  return tag<GlobalTag>(std::move(fullVector));
}

template <class Policy>
global<Vector<Policy>> Vector<Policy>::fromDistributedInCanonicalOrder(
    const distributed<Vector> &vector, const Mesh<Policy> &mesh,
    const size_type vertexDof, const size_type elementDof) {
  const auto localNumberOfVertices = mesh.localNumberOfVertices();
  const auto totalNumberOfVertices = mesh.totalNumberOfVertices();
  const auto localNumberOfElements = mesh.localNumberOfElements();

  std::vector<size_type> fromIndices;
  fromIndices.reserve(vertexDof * localNumberOfVertices +
                      elementDof * localNumberOfElements);
  std::vector<size_type> toIndices;
  toIndices.reserve(fromIndices.size());

  for (const auto &vertex : mesh.localVertices()) {
    const auto range = vertex.globalDofLineRange();
    assert(std::abs(range.second - range.first) == vertexDof &&
           "The number of dofs must be the same for all vertices.");

    for (auto i = range.first; i < range.second; ++i) {
      fromIndices.push_back(i);
      toIndices.push_back(vertex.index() * vertexDof + (i - range.first));
    }
  }
  for (const auto &element : mesh.localElements()) {
    const auto range = element.globalDofLineRange();
    assert(range.second - range.first == elementDof &&
           "The number of dofs must be the same for all elements.");

    for (auto i = range.first; i < range.second; ++i) {
      fromIndices.push_back(i);
      toIndices.push_back(vertexDof * totalNumberOfVertices +
                          element.index() * elementDof + (i - range.first));
    }
  }

  const auto toGlobalIS = [](const std::vector<size_type> &indices) {
    const auto isLocal = [&indices]() {
      auto is = IS();
      Policy::handleError(ISCreateGeneral(Policy::communicator(),
                                          indices.size(), indices.data(),
                                          PETSC_USE_POINTER, &is));
      return makeUniqueEntity<Policy>(is);
    }();

    auto isGlobal = IS();
    ISAllGather(isLocal.get(), &isGlobal);
    return makeUniqueEntity<Policy>(isGlobal);
  };

  const auto globalFrom = toGlobalIS(fromIndices);
  const auto globalTo = toGlobalIS(toIndices);

  auto fullVector = fromDistributed(vector);

  Policy::handleError(
      VecPermute(fullVector.unwrap().data(), globalFrom.get(), PETSC_FALSE));
  Policy::handleError(
      VecPermute(fullVector.unwrap().data(), globalTo.get(), PETSC_TRUE));

  return fullVector;
}

template <class Policy>
Vector<Policy> Vector<Policy>::fromLayoutOf(const Vector &vector) {
  auto vec = Vec{};
  Policy::handleError(VecDuplicate(vector.data(), &vec));
  auto newVector = Vector(makeUniqueEntity<Policy>(vec));
  newVector.setZero();
  return newVector;
}

template <class Policy>
template <InsertMode Mode>
class Vector<Policy>::Inserter {
public:
  using value_type = Vector<Policy>::value_type;
  using size_type = Vector<Policy>::size_type;

  Inserter(Vector *const vector) : _vector(vector) {}

  size_type size() const { return _vector->size(); }

  const Inserter &elements(const std::vector<size_type> &indices,
                           const std::vector<value_type> &values) const {
    Policy::handleError(VecSetValues(_vector->_vec.get(), indices.size(),
                                     indices.data(), values.data(), Mode));
    return *this;
  }

  InsertionProxy<Mode>
  operator()(const typename Vector::size_type index) const {
    return InsertionProxy<Mode>{
        [this, index](const value_type value) { this->element(index, value); }};
  }

  const Inserter &element(const size_type index, const value_type value) const {
    Policy::handleError(VecSetValue(_vector->_vec.get(), index, value, Mode));
    return *this;
  }

  void fill(const value_type value) const {
    flush();
    _vector->fill(value);
  }

  void setZero() const {
    flush();
    _vector->setZero();
  }

  void flush() const {
    Policy::handleError(VecAssemblyBegin(_vector->_vec.get()));
    Policy::handleError(VecAssemblyEnd(_vector->_vec.get()));
  }

  ~Inserter() { flush(); }

private:
  Vector *_vector = nullptr;
};

template <class Policy>
typename Vector<Policy>::template Inserter<INSERT_VALUES>
Vector<Policy>::replace() {
  return Inserter<INSERT_VALUES>(this);
}

template <class Policy>
typename Vector<Policy>::template Inserter<ADD_VALUES> Vector<Policy>::add() {
  return Inserter<ADD_VALUES>(this);
}

template <class Policy>
Vector<Policy>
Vector<Policy>::fromList(const std::initializer_list<value_type> list) {
  auto vector = Vector(list.size());

  const std::vector<value_type> values(list.begin(), list.end());

  std::vector<size_type> indices(list.size());
  std::iota(indices.begin(), indices.end(), size_type{0});

  vector.replace().elements(indices, values);
  return vector;
}

template <class Policy>
typename Vector<Policy>::value_type Vector<Policy>::norm() const {
  auto result = value_type{};
  Policy::handleError(VecNorm(_vec.get(), NORM_2, &result));
  return result;
}

template <class Policy> void Vector<Policy>::print() const {
  Policy::handleError(VecView(_vec.get(), nullptr));
}

template <class Policy>
typename Vector<Policy>::value_type Vector<Policy>::
operator()(const size_type index) const {
  auto value = value_type{0.};
  size_type indices[] = {index};
  Policy::handleError(VecGetValues(_vec.get(), 1, indices, &value));
  return value;
}

template <class Policy>
typename Vector<Policy>::value_type Vector<Policy>::
operator[](const size_type index) const {
  return (*this)(index);
}

template <class Policy> void Vector<Policy>::scale(const value_type factor) {
  Policy::handleError(VecScale(_vec.get(), factor));
}

template <class Policy>
void Vector<Policy>::timesAlphaPlusBetaX(const value_type alpha,
                                         const value_type beta,
                                         const Vector &x) {
  Policy::handleError(VecAXPBY(_vec.get(), beta, alpha, x.data()));
}

template <class Policy>
void Vector<Policy>::timesAlphaPlusBetaXPlusGammaY(const value_type alpha,
                                                   const value_type beta,
                                                   const Vector &x,
                                                   const value_type gamma,
                                                   const Vector &y) {
  Policy::handleError(
      VecAXPBYPCZ(_vec.get(), beta, gamma, alpha, x.data(), y.data()));
}

template <class Policy>
typename Vector<Policy>::size_type Vector<Policy>::size() const {
  size_type size = 0;
  Policy::handleError(VecGetSize(_vec.get(), &size));
  return size;
}

template <class Policy> void Vector<Policy>::setZero() {
  Policy::handleError(VecZeroEntries(_vec.get()));
}

template <class Policy>
std::pair<typename Vector<Policy>::size_type,
          typename Vector<Policy>::size_type>
Vector<Policy>::localRowRange() const {
  auto return_value = std::pair<size_type, size_type>();
  Policy::handleError(VecGetOwnershipRange(_vec.get(), &return_value.first,
                                           &return_value.second));
  return return_value;
}

template <class Policy>
typename Vector<Policy>::const_iterator Vector<Policy>::localBegin() const {
  return const_iterator{this, localRowRange().first};
}

template <class Policy>
typename Vector<Policy>::const_iterator Vector<Policy>::localEnd() const {
  return const_iterator{this, localRowRange().second};
}

template <class Policy>
IteratorRange<typename Vector<Policy>::const_iterator>
Vector<Policy>::localRange() const {
  const auto range = localRowRange();
  return {const_iterator{this, range.first},
          const_iterator{this, range.second}};
}

template <class Policy> void Vector<Policy>::fill(const value_type value) {
  Policy::handleError(VecSet(_vec.get(), value));
}

template <class Policy> Vec Vector<Policy>::data() const { return _vec.get(); }

template <class Policy>
Vector<Policy> Vector<Policy>::clone(const Vector &vector) {
  auto result = fromLayoutOf(vector);
  Policy::handleError(VecCopy(vector.data(), result.data()));
  return result;
}

template <class Policy>
void Vector<Policy>::addAx(const matrix_type &A, const distributed<Vector> &x) {
  Policy::handleError(
      MatMultAdd(A.data(), x.unwrap().data(), _vec.get(), _vec.get()));
}
} // namespace cpppetsc
} // namespace ae108