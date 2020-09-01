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

#include "ae108/cpppetsc/TaggedVector.h"
#include <utility>
#include <vector>

namespace ae108 {
namespace cpppetsc {

template <class IteratorType> class LocalVertexView {
  friend IteratorType;

public:
  using mesh_type = typename IteratorType::mesh_type;

  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;
  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;

  explicit LocalVertexView();

  /**
   * @param mesh A valid pointer to a mesh_type instance.
   * @param vertexPointIndex A valid vertex point index.
   */
  explicit LocalVertexView(const mesh_type *const mesh,
                           const size_type vertexPointIndex);

  /**
   * @brief Converts the vertex view to a vertex view for a cloned mesh.
   * @param mesh is a Mesh that was cloned from the Mesh that this view refers
   * to.
   */
  LocalVertexView forClonedMesh(const mesh_type *const mesh) const;

  /**
   * @brief Returns the global vertex index.
   */
  size_type index() const;

  /**
   * @brief Returns the number of degrees of freedom associated with the vertex.
   */
  size_type numberOfDofs() const;

  /**
   * @brief Returns the local range of lines (e.g. in a vector) that is
   * associated with the degrees of freedom of this vertex in the following
   * format: (first line, last line + 1).
   */
  std::pair<size_type, size_type> localDofLineRange() const;

  /**
   * @brief Returns the global range of lines (e.g. in a vector) that is
   * associated with the degrees of freedom of this vertex in the following
   * format: (first line, last line + 1).
   *
   * @remark Returns negative line numbers (-begin - 1, -end - 1) if the point
   * data is not owned locally (e.g. for ghost points).
   */
  std::pair<size_type, size_type> globalDofLineRange() const;

  /**
   * @brief Copies the data belonging to the vertex to data.
   */
  void copyVertexData(const local<vector_type> &vector,
                      std::vector<value_type> *data) const;

  /**
   * @brief Adds the provided vertex data to the vector.
   */
  void addVertexData(const std::vector<value_type> &data,
                     local<vector_type> *const vector) const;

  /**
   * @brief Sets the provided vertex data in the vector.
   */
  void setVertexData(const std::vector<value_type> &data,
                     local<vector_type> *const vector) const;

  /**
   * @brief Adds the provided matrix data to the matrix.
   */
  void addVertexMatrix(const std::vector<value_type> &data,
                       matrix_type *const matrix) const;

  friend bool operator==(const LocalVertexView &lhs,
                         const LocalVertexView &rhs) {
    return lhs._mesh == rhs._mesh && lhs._entityIndex == rhs._entityIndex;
  }

private:
  const mesh_type *_mesh = nullptr;
  size_type _entityIndex = 0;
};
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cpppetsc {

template <class IteratorType>
LocalVertexView<IteratorType>::LocalVertexView() = default;

template <class IteratorType>
LocalVertexView<IteratorType>::LocalVertexView(const mesh_type *const mesh,
                                               const size_type vertexPointIndex)
    : _mesh(mesh), _entityIndex(vertexPointIndex) {}

template <class IteratorType>
LocalVertexView<IteratorType> LocalVertexView<IteratorType>::forClonedMesh(
    const mesh_type *const mesh) const {
  return LocalVertexView{mesh, _entityIndex};
}

template <class IteratorType>
typename LocalVertexView<IteratorType>::size_type
LocalVertexView<IteratorType>::index() const {
  return _mesh->vertexPointIndexToGlobalIndex(_entityIndex);
}

template <class IteratorType>
typename LocalVertexView<IteratorType>::size_type
LocalVertexView<IteratorType>::numberOfDofs() const {
  return _mesh->numberOfDofs(_entityIndex);
}
template <class IteratorType>
std::pair<typename LocalVertexView<IteratorType>::size_type,
          typename LocalVertexView<IteratorType>::size_type>
LocalVertexView<IteratorType>::localDofLineRange() const {
  return _mesh->localDofLineRange(_entityIndex);
}

template <class IteratorType>
std::pair<typename LocalVertexView<IteratorType>::size_type,
          typename LocalVertexView<IteratorType>::size_type>
LocalVertexView<IteratorType>::globalDofLineRange() const {
  return _mesh->globalDofLineRange(_entityIndex);
}

template <class IteratorType>
void LocalVertexView<IteratorType>::copyVertexData(
    const local<vector_type> &vector, std::vector<value_type> *data) const {
  return _mesh->copyEntityData(_entityIndex, vector, data);
}

template <class IteratorType>
void LocalVertexView<IteratorType>::addVertexData(
    const std::vector<value_type> &data,
    local<vector_type> *const vector) const {
  return _mesh->addEntityData(_entityIndex, data, vector);
}

template <class IteratorType>
void LocalVertexView<IteratorType>::setVertexData(
    const std::vector<value_type> &data,
    local<vector_type> *const vector) const {
  return _mesh->setEntityData(_entityIndex, data, vector);
}

template <class IteratorType>
void LocalVertexView<IteratorType>::addVertexMatrix(
    const std::vector<value_type> &data, matrix_type *const matrix) const {
  return _mesh->addEntityMatrix(_entityIndex, data, matrix);
}
} // namespace cpppetsc
} // namespace ae108