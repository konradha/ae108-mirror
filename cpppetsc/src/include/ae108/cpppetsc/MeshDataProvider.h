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

#include "ae108/cpppetsc/SharedEntity.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include <petscdm.h>
#include <utility>
#include <vector>

namespace ae108 {
namespace cpppetsc {

/**
 * brief Provides the necessary data for local entity views.
 */
template <class MeshType_> class MeshDataProvider {
public:
  using mesh_type = MeshType_;

  using policy_type = typename mesh_type::policy_type;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;
  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;

  explicit MeshDataProvider(SharedEntity<DM, policy_type> mesh,
                            const size_type totalNumberOfElements);

  /**
   * @brief Convert an internal element index to a 'canonical' element index.
   */
  size_type elementPointIndexToGlobalIndex(const size_type pointIndex) const;

  /**
   * @brief Convert an internal vertex index to a 'canonical' vertex index.
   */
  size_type vertexPointIndexToGlobalIndex(const size_type pointIndex) const;

  /**
   * @brief Returns the number of vertices of the element.
   */
  size_type numberOfVertices(const size_type elementPointIndex) const;

  /**
   * @brief Returns the vertex indices for the element.
   */
  std::vector<size_type> vertices(const size_type elementPointIndex) const;

  /**
   * @brief Returns the global line range that corresponds to the degrees of
   * freedom of the vertex or element.
   */
  std::pair<size_type, size_type>
  localDofLineRange(const size_type entityPointIndex) const;

  /**
   * @brief Returns the global line range that corresponds to the degrees of
   * freedom of the vertex or element.
   *
   * @remark Returns negative line numbers (-begin - 1, -end - 1) if the point
   * data is not owned locally (e.g. for ghost points). In particular, first >=
   * end in this case.
   */
  std::pair<size_type, size_type>
  globalDofLineRange(const size_type entityPointIndex) const;

  size_type numberOfDofs(const size_type entityPointIndex) const;

  /**
   * @brief Copies data corresponding to the specified entity from the vector
   * to the data buffer.
   *
   * @param entityPointIndex Specifies the entity.
   * @param localVector A mesh-local vector.
   * @param data A valid pointer to a buffer object.
   */
  void copyEntityData(const size_type entityPointIndex,
                      const local<vector_type> &localVector,
                      std::vector<value_type> *const data) const;

  /**
   * @brief Add data corresponding to the specified entity from the data buffer
   * to the vector.
   *
   * @param entityPointIndex Specifies the entity.
   * @param data The data buffer.
   * @param localVector A valid pointer to a mesh-local vector.
   */
  void addEntityData(const size_type entityPointIndex,
                     const std::vector<value_type> &data,
                     local<vector_type> *const localVector) const;

  /**
   * @brief Set data corresponding to the specified entity from the data buffer
   * in the vector.
   *
   * @param entityPointIndex Specifies the entity.
   * @param data The data buffer.
   * @param localVector A valid pointer to a mesh-local vector.
   */
  void setEntityData(const size_type entityPointIndex,
                     const std::vector<value_type> &data,
                     local<vector_type> *const localVector) const;

  /**
   * @brief Add data corresponding to the specified entity from the data buffer
   * to the matrix.
   *
   * @param entityPointIndex Specifies the entity.
   * @param data The data buffer.
   * @param matrix A valid pointer to a matrix.
   */
  void addEntityMatrix(const size_type entityPointIndex,
                       const std::vector<value_type> &data,
                       matrix_type *const matrix) const;

  /**
   * @brief Asserts that the vector has been created with the associated mesh.
   */
  void assertCorrectBaseMesh(const vector_type &vector) const;

private:
  SharedEntity<DM, policy_type> _mesh;
  size_type _totalNumberOfElements;
};

/**
 * @brief Creates a data provider from the given mesh.
 * @param mesh Valid nonzero pointer.
 */
template <class MeshType>
MeshDataProvider<MeshType>
createDataProviderFromMesh(const MeshType *const mesh);

} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include <algorithm>
#include <cassert>
#include <petscdmplex.h>
#include <petscsf.h>

namespace ae108 {
namespace cpppetsc {

template <class MeshType_>
MeshDataProvider<MeshType_>::MeshDataProvider(
    SharedEntity<DM, policy_type> mesh, const size_type totalNumberOfElements)
    : _mesh(std::move(mesh)), _totalNumberOfElements(totalNumberOfElements) {}

template <class MeshType_>
typename MeshDataProvider<MeshType_>::size_type
MeshDataProvider<MeshType_>::elementPointIndexToGlobalIndex(
    const size_type pointIndex) const {
  auto migration = PetscSF();
  policy_type::handleError(DMPlexGetMigrationSF(_mesh.get(), &migration));
  if (migration) {
    const PetscSFNode *nodes = nullptr;
    policy_type::handleError(
        PetscSFGetGraph(migration, nullptr, nullptr, nullptr, &nodes));
    return nodes[pointIndex].index;
  }

  return pointIndex;
}

template <class MeshType_>
typename MeshDataProvider<MeshType_>::size_type
MeshDataProvider<MeshType_>::vertexPointIndexToGlobalIndex(
    const size_type pointIndex) const {
  return elementPointIndexToGlobalIndex(pointIndex) - _totalNumberOfElements;
}

template <class MeshType_>
typename MeshDataProvider<MeshType_>::size_type
MeshDataProvider<MeshType_>::numberOfVertices(
    const size_type elementPointIndex) const {
  size_type size = 0;
  policy_type::handleError(
      DMPlexGetConeSize(_mesh.get(), elementPointIndex, &size));
  return size;
}

template <class MeshType_>
std::vector<typename MeshDataProvider<MeshType_>::size_type>
MeshDataProvider<MeshType_>::vertices(const size_type elementPointIndex) const {
  const auto size = numberOfVertices(elementPointIndex);
  const size_type *cone = nullptr;
  policy_type::handleError(
      DMPlexGetCone(_mesh.get(), elementPointIndex, &cone));
  std::vector<size_type> output(size);
  std::transform(cone, cone + size, output.begin(),
                 [this](const size_type entityPointIndex) {
                   return vertexPointIndexToGlobalIndex(entityPointIndex);
                 });
  return output;
}

template <class MeshType_>
std::pair<typename MeshDataProvider<MeshType_>::size_type,
          typename MeshDataProvider<MeshType_>::size_type>
MeshDataProvider<MeshType_>::localDofLineRange(
    const size_type entityPointIndex) const {
  auto output = std::pair<size_type, size_type>();
  policy_type::handleError(DMPlexGetPointLocal(_mesh.get(), entityPointIndex,
                                               &output.first, &output.second));
  return output;
}

template <class MeshType_>
std::pair<typename MeshDataProvider<MeshType_>::size_type,
          typename MeshDataProvider<MeshType_>::size_type>
MeshDataProvider<MeshType_>::globalDofLineRange(
    const size_type entityPointIndex) const {
  auto output = std::pair<size_type, size_type>();
  policy_type::handleError(DMPlexGetPointGlobal(_mesh.get(), entityPointIndex,
                                                &output.first, &output.second));
  return output;
}

template <class MeshType_>
typename MeshDataProvider<MeshType_>::size_type
MeshDataProvider<MeshType_>::numberOfDofs(
    const size_type entityPointIndex) const {
  auto output = size_type{0};
  auto section = PetscSection();
  policy_type::handleError(DMGetSection(_mesh.get(), &section));

  policy_type::handleError(
      PetscSectionGetDof(section, entityPointIndex, &output));

  return output;
}

template <class MeshType_>
void MeshDataProvider<MeshType_>::copyEntityData(
    const size_type entityPointIndex, const local<vector_type> &localVector,
    std::vector<value_type> *const data) const {
  assertCorrectBaseMesh(localVector.unwrap());
  auto size = size_type{0};
  policy_type::handleError(
      DMPlexVecGetClosure(_mesh.get(), nullptr, localVector.unwrap().data(),
                          entityPointIndex, &size, nullptr));
  data->resize(size);
  if (size == 0)
    return;

  auto dataPtr = data->data();
  policy_type::handleError(
      DMPlexVecGetClosure(_mesh.get(), nullptr, localVector.unwrap().data(),
                          entityPointIndex, &size, &dataPtr));
}

template <class MeshType_>
void MeshDataProvider<MeshType_>::addEntityData(
    const size_type entityPointIndex, const std::vector<value_type> &data,
    local<vector_type> *const localVector) const {
  assertCorrectBaseMesh(localVector->unwrap());
  policy_type::handleError(
      DMPlexVecSetClosure(_mesh.get(), nullptr, localVector->unwrap().data(),
                          entityPointIndex, data.data(), ADD_VALUES));
}

template <class MeshType_>
void MeshDataProvider<MeshType_>::setEntityData(
    const size_type entityPointIndex, const std::vector<value_type> &data,
    local<vector_type> *const localVector) const {
  assertCorrectBaseMesh(localVector->unwrap());
  policy_type::handleError(
      DMPlexVecSetClosure(_mesh.get(), nullptr, localVector->unwrap().data(),
                          entityPointIndex, data.data(), INSERT_VALUES));
}

template <class MeshType_>
void MeshDataProvider<MeshType_>::addEntityMatrix(
    const size_type entityPointIndex, const std::vector<value_type> &data,
    matrix_type *const matrix) const {
  policy_type::handleError(DMPlexMatSetClosure(_mesh.get(), nullptr, nullptr,
                                               matrix->data(), entityPointIndex,
                                               data.data(), ADD_VALUES));
}

template <class MeshType_>
void MeshDataProvider<MeshType_>::assertCorrectBaseMesh(
    const vector_type &vector) const {
#ifndef NDEBUG
  const auto extractDM = [](const vector_type &vector) {
    auto dm = DM{};
    policy_type::handleError(VecGetDM(vector.data(), &dm));
    return dm;
  };
  const auto dm = extractDM(vector);
  assert((!dm || dm == _mesh.get()) &&
         "The vector was created using a different mesh.");
#else
  static_cast<void>(vector);
#endif
}

template <class MeshType>
MeshDataProvider<MeshType>
createDataProviderFromMesh(const MeshType *const mesh) {
  assert(mesh);
  return MeshDataProvider<MeshType>{
      makeSharedEntity<typename MeshType::policy_type>(mesh->data(),
                                                       OwnershipType::Share),
      mesh->totalNumberOfElements(),
  };
}
} // namespace cpppetsc
} // namespace ae108