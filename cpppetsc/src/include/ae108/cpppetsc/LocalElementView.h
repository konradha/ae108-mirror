// © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#include "ae108/cpppetsc/MeshDataProvider.h"
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

  using DataProvider = MeshDataProvider<mesh_type>;

  /**
   * @param data Provides mesh information to the view.
   * @param elementPointIndex A valid element point index.
   */
  explicit LocalElementView(DataProvider data,
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

private:
  DataProvider _data;
  size_type _entityIndex;
};
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include <cassert>

namespace ae108 {
namespace cpppetsc {

template <class IteratorType>
LocalElementView<IteratorType>::LocalElementView(
    DataProvider data, const size_type elementPointIndex)
    : _data(std::move(data)), _entityIndex(elementPointIndex) {}

template <class IteratorType>
LocalElementView<IteratorType> LocalElementView<IteratorType>::forClonedMesh(
    const mesh_type *const mesh) const {
  assert(mesh);
  return LocalElementView{createDataProviderFromMesh(mesh), _entityIndex};
}

template <class IteratorType>
std::vector<typename LocalElementView<IteratorType>::size_type>
LocalElementView<IteratorType>::vertexIndices() const {
  return _data.vertices(_entityIndex);
}

template <class IteratorType>
typename LocalElementView<IteratorType>::size_type
LocalElementView<IteratorType>::index() const {
  return _data.elementPointIndexToGlobalIndex(_entityIndex);
}

template <class IteratorType>
typename LocalElementView<IteratorType>::size_type
LocalElementView<IteratorType>::numberOfVertices() const {
  return _data.numberOfVertices(_entityIndex);
}

template <class IteratorType>
typename LocalElementView<IteratorType>::size_type
LocalElementView<IteratorType>::numberOfDofs() const {
  return _data.numberOfDofs(_entityIndex);
}

template <class IteratorType>
std::pair<typename LocalElementView<IteratorType>::size_type,
          typename LocalElementView<IteratorType>::size_type>
LocalElementView<IteratorType>::localDofLineRange() const {
  return _data.localDofLineRange(_entityIndex);
}

template <class IteratorType>
std::pair<typename LocalElementView<IteratorType>::size_type,
          typename LocalElementView<IteratorType>::size_type>
LocalElementView<IteratorType>::globalDofLineRange() const {
  return _data.globalDofLineRange(_entityIndex);
}

template <class IteratorType>
void LocalElementView<IteratorType>::copyElementData(
    const local<vector_type> &vector, std::vector<value_type> *data) const {
  _data.copyEntityData(_entityIndex, vector, data);
}

template <class IteratorType>
void LocalElementView<IteratorType>::addElementData(
    const std::vector<value_type> &data,
    local<vector_type> *const vector) const {
  _data.addEntityData(_entityIndex, data, vector);
}

template <class IteratorType>
void LocalElementView<IteratorType>::setElementData(
    const std::vector<value_type> &data,
    local<vector_type> *const vector) const {
  _data.setEntityData(_entityIndex, data, vector);
}

template <class IteratorType>
void LocalElementView<IteratorType>::addElementMatrix(
    const std::vector<value_type> &data, matrix_type *const matrix) const {
  _data.addEntityMatrix(_entityIndex, data, matrix);
}
} // namespace cpppetsc
} // namespace ae108