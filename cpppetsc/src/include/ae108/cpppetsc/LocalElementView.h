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
#include <vector>

namespace ae108 {
namespace cpppetsc {

template <class MeshType_> class LocalElementView {
public:
  using mesh_type = MeshType_;

  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;
  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;

  explicit LocalElementView();

  /**
   * @param mesh A valid pointer to a mesh_type instance.
   * @param elementPointIndex A valid element point index.
   */
  explicit LocalElementView(const mesh_type *const mesh,
                            const size_type elementPointIndex);

  /**
   * @brief Converts the element view to an element view for a cloned mesh.
   * @param mesh is a Mesh that was cloned from the Mesh that this view refers
   * to.
   */
  LocalElementView forClonedMesh(const mesh_type *const mesh) const;

  /**
   * @brief Returns the global element index.
   */
  size_type index() const;

  /**
   * @brief Returns the number of vertices of the element.
   */
  size_type numberOfVertices() const;

  /**
   * @brief Returns the number of degrees of freedom of the element. This number
   * does not include the degrees of freedom of the associated vertices.
   */
  size_type numberOfDofs() const;

  /**
   * @brief Returns the local range of lines (e.g. in a vector) that is
   * associated with the degrees of freedom of this element in the following
   * format: (first line, last line + 1).
   */
  std::pair<size_type, size_type> localDofLineRange() const;

  /**
   * @brief Returns the global range of lines (e.g. in a vector) that is
   * associated with the degrees of freedom of this element in the following
   * format: (first line, last line + 1).
   */
  std::pair<size_type, size_type> globalDofLineRange() const;

  /**
   * @brief The global vertex indices of the element.
   */
  std::vector<size_type> vertexIndices() const;

  /**
   * @brief Copies the data belonging to the element to data.
   */
  void copyElementData(const local<vector_type> &vector,
                       std::vector<value_type> *data) const;

  /**
   * @brief Adds the provided element data to the vector.
   */
  void addElementData(const std::vector<value_type> &data,
                      local<vector_type> *const vector) const;

  /**
   * @brief Sets the provided element data in the vector.
   */
  void setElementData(const std::vector<value_type> &data,
                      local<vector_type> *const vector) const;

  /**
   * @brief Adds the provided matrix data to the matrix.
   */
  void addElementMatrix(const std::vector<value_type> &data,
                        matrix_type *const matrix) const;

  friend bool operator==(const LocalElementView &lhs,
                         const LocalElementView &rhs) {
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
LocalElementView<IteratorType>::LocalElementView() = default;

template <class IteratorType>
LocalElementView<IteratorType>::LocalElementView(
    const mesh_type *const mesh, const size_type elementPointIndex)
    : _mesh(mesh), _entityIndex(elementPointIndex) {}

template <class IteratorType>
LocalElementView<IteratorType> LocalElementView<IteratorType>::forClonedMesh(
    const mesh_type *const mesh) const {
  return LocalElementView{mesh, _entityIndex};
}

template <class IteratorType>
std::vector<typename LocalElementView<IteratorType>::size_type>
LocalElementView<IteratorType>::vertexIndices() const {
  return _mesh->vertices(_entityIndex);
}

template <class IteratorType>
typename LocalElementView<IteratorType>::size_type
LocalElementView<IteratorType>::index() const {
  return _mesh->elementPointIndexToGlobalIndex(_entityIndex);
}

template <class IteratorType>
typename LocalElementView<IteratorType>::size_type
LocalElementView<IteratorType>::numberOfVertices() const {
  return _mesh->numberOfVertices(_entityIndex);
}

template <class IteratorType>
typename LocalElementView<IteratorType>::size_type
LocalElementView<IteratorType>::numberOfDofs() const {
  return _mesh->numberOfDofs(_entityIndex);
}

template <class IteratorType>
std::pair<typename LocalElementView<IteratorType>::size_type,
          typename LocalElementView<IteratorType>::size_type>
LocalElementView<IteratorType>::localDofLineRange() const {
  return _mesh->localDofLineRange(_entityIndex);
}

template <class IteratorType>
std::pair<typename LocalElementView<IteratorType>::size_type,
          typename LocalElementView<IteratorType>::size_type>
LocalElementView<IteratorType>::globalDofLineRange() const {
  return _mesh->globalDofLineRange(_entityIndex);
}

template <class IteratorType>
void LocalElementView<IteratorType>::copyElementData(
    const local<vector_type> &vector, std::vector<value_type> *data) const {
  _mesh->copyEntityData(_entityIndex, vector, data);
}

template <class IteratorType>
void LocalElementView<IteratorType>::addElementData(
    const std::vector<value_type> &data,
    local<vector_type> *const vector) const {
  _mesh->addEntityData(_entityIndex, data, vector);
}

template <class IteratorType>
void LocalElementView<IteratorType>::setElementData(
    const std::vector<value_type> &data,
    local<vector_type> *const vector) const {
  _mesh->setEntityData(_entityIndex, data, vector);
}

template <class IteratorType>
void LocalElementView<IteratorType>::addElementMatrix(
    const std::vector<value_type> &data, matrix_type *const matrix) const {
  _mesh->addEntityMatrix(_entityIndex, data, matrix);
}
} // namespace cpppetsc
} // namespace ae108