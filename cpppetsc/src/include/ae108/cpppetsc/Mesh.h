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

#include "ae108/cpppetsc/IteratorRange.h"
#include "ae108/cpppetsc/LocalElementIterator.h"
#include "ae108/cpppetsc/LocalElementView.h"
#include "ae108/cpppetsc/LocalVertexIterator.h"
#include "ae108/cpppetsc/LocalVertexView.h"
#include "ae108/cpppetsc/Matrix_fwd.h"
#include "ae108/cpppetsc/Mesh_fwd.h"
#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/UniqueEntity.h"
#include "ae108/cpppetsc/Vector_fwd.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscis.h>
#include <petscmath.h>
#include <petscsf.h>
#include <petscsys.h>
#include <utility>
#include <vector>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class Mesh {
  friend LocalElementView<LocalElementIterator<Mesh>>;
  friend LocalVertexView<LocalVertexIterator<Mesh>>;

public:
  using size_type = PetscInt;
  using value_type = PetscScalar;
  using vector_type = Vector<Policy>;
  using matrix_type = Matrix<Policy>;

  static constexpr size_type IGNORE_VERTEX_INDEX = -1;

  /**
   * @brief Calls the second overload with topological and coordinate
   * dimension both equal to `dimension`.
   */
  template <class Container>
  static Mesh fromConnectivity(const size_type dimension,
                               const Container &elementVertexIDs,
                               const size_type numberOfVertices,
                               const size_type dofPerVertex,
                               const size_type dofPerElement = 0);

  using TopologicalDimension =
      TaggedEntity<size_type, struct TopologicalDimensionTag>;
  using CoordinateDimension =
      TaggedEntity<size_type, struct CoordinateDimensionTag>;

  /**
   * @brief Creates a Mesh from the connectivity representation.
   *
   * @tparam Container Requirements: size(); Container::value_type requirements:
   * size(), data().
   *
   * @param elementVertexIDs A container of containers: the container at index i
   * contains the vertex IDs of the element with index i. Each vertex ID must
   * lie in [0, numberOfVertices) u {IGNORE_VERTEX_INDEX}. IDs of value
   * IGNORE_VERTEX_INDEX are ignored.
   *
   * @param numberOfVertices The maximum vertex ID in the mesh.
   */
  template <class Container>
  static Mesh fromConnectivity(const TopologicalDimension topologicalDimension,
                               const CoordinateDimension coordinateDimension,
                               const Container &elementVertexIDs,
                               const size_type numberOfVertices,
                               const size_type dofPerVertex,
                               const size_type dofPerElement);

  using const_element_iterator = LocalElementIterator<Mesh>;
  using const_iterator = const_element_iterator;

  /**
   * @brief Clone the mesh with a different default section.
   */
  Mesh cloneWithDofs(const size_type dofPerVertex,
                     const size_type dofPerElement) const;

  /**
   * @brief Makes it possible to iterate over (views of) the local elements.
   */
  IteratorRange<const_iterator> localElements() const;

  using const_vertex_iterator = LocalVertexIterator<Mesh>;
  /**
   * @brief Makes it possible to iterate over (views of) the local vertices.
   */
  IteratorRange<const_vertex_iterator> localVertices() const;

  /**
   * @brief Returns the total number of elements in the mesh.
   */
  size_type totalNumberOfElements() const;

  /**
   * @brief Returns the local number of elements in the mesh.
   */
  size_type localNumberOfElements() const;

  /**
   * @brief Returns the total number of vertices in the mesh.
   */
  size_type totalNumberOfVertices() const;

  /**
   * @brief Returns the local number of vertices in the mesh.
   */
  size_type localNumberOfVertices() const;

  /**
   * @brief Returns the dimension of the elements.
   */
  TopologicalDimension topologicalDimension() const;

  /**
   * @brief Returns the dimension of the coordinates.
   */
  CoordinateDimension coordinateDimension() const;

  /**
   * @brief Fill the local vector with the corresponding values from the
   * global vector.
   */
  void copyToLocalVector(const distributed<vector_type> &globalVector,
                         local<vector_type> *const localVector) const;

  /**
   * @brief Add the values in the local vector to the corresponding values in
   * the global vector.
   */
  void addToGlobalVector(const local<vector_type> &localVector,
                         distributed<vector_type> *const globalVector) const;

  /**
   * @brief Copy the values in the local vector to the corresponding values
   * in the global vector.
   *
   * @remark Be careful with overlapping areas (e.g. for ghost points) on
   * different nodes.
   */
  void copyToGlobalVector(const local<vector_type> &localVector,
                          distributed<vector_type> *const globalVector) const;

  /**
   * @brief Prints the mesh to stdout.
   */
  void print() const;

  /**
   * @brief Returns the internal mesh.
   *
   * @remark For internal use only.
   */
  DM data() const;

private:
  /**
   * @brief Creates a DM and sets the provided members.
   */
  explicit Mesh(const size_type totalNumberOfElements,
                const size_type totalNumberOfVertices);

  static UniqueEntity<DM> createDM();

  /**
   * @brief Adds the default uniform section to the mesh.
   */
  void addSection(const size_type dofPerVertex, const size_type dofPerElement);

  /**
   * @brief Adds the chart as described by the connectivity to the mesh.
   */
  template <class Container>
  static void addChart(const Container &elementVertexIDs,
                       const size_type numberOfVertices, const DM dm);

  /**
   * @brief Sets the adjacency rules of the provided mesh.
   */
  static void setAdjacencyRules(const DM dm);

  /**
   * @brief Distributes the provided mesh over the MPI nodes.
   */
  static void distributeMesh(Mesh *const);

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
   * @brief Asserts that the vector has been created with this mesh.
   */
  void assertCorrectBaseMesh(const vector_type &vector) const;

  UniqueEntity<DM> _mesh;
  size_type _totalNumberOfElements = 0;
  size_type _totalNumberOfVertices = 0;
};

template <class Policy>
constexpr typename Mesh<Policy>::size_type Mesh<Policy>::IGNORE_VERTEX_INDEX;

extern template class Mesh<SequentialComputePolicy>;
extern template class Mesh<ParallelComputePolicy>;
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

template <class Policy>
IteratorRange<typename Mesh<Policy>::const_iterator>
Mesh<Policy>::localElements() const {
  auto start = size_type{0};
  auto stop = size_type{0};
  Policy::handleError(DMPlexGetHeightStratum(_mesh.get(), 0, &start, &stop));
  return IteratorRange<const_iterator>{const_iterator{this, start},
                                       const_iterator{this, stop}};
}

template <class Policy>
IteratorRange<typename Mesh<Policy>::const_vertex_iterator>
Mesh<Policy>::localVertices() const {
  auto start = size_type{0};
  auto stop = size_type{0};
  Policy::handleError(DMPlexGetDepthStratum(_mesh.get(), 0, &start, &stop));
  return IteratorRange<const_vertex_iterator>{
      const_vertex_iterator{this, start}, const_vertex_iterator{this, stop}};
}

template <class Policy>
template <class Container>
Mesh<Policy> Mesh<Policy>::fromConnectivity(const size_type dimension,
                                            const Container &elementVertexIDs,
                                            const size_type numberOfVertices,
                                            const size_type dofPerVertex,
                                            const size_type dofPerElement) {
  return fromConnectivity(TopologicalDimension(dimension),
                          CoordinateDimension(dimension), elementVertexIDs,
                          numberOfVertices, dofPerVertex, dofPerElement);
}

template <class Policy>
template <class Container>
Mesh<Policy> Mesh<Policy>::fromConnectivity(
    const TopologicalDimension topologicalDimension,
    const CoordinateDimension coordinateDimension,
    const Container &elementVertexIDs, const size_type numberOfVertices,
    const size_type dofPerVertex, const size_type dofPerElement) {
  auto mesh =
      Mesh(static_cast<size_type>(elementVertexIDs.size()), numberOfVertices);

  Policy::handleError(
      DMSetCoordinateDim(mesh._mesh.get(), coordinateDimension));
  Policy::handleError(DMSetDimension(mesh._mesh.get(), topologicalDimension));

  addChart(elementVertexIDs, numberOfVertices, mesh._mesh.get());

  // calculate the strata
  Policy::handleError(DMPlexStratify(mesh._mesh.get()));

  distributeMesh(&mesh);

  mesh.addSection(dofPerVertex, dofPerElement);

  return mesh;
}

template <class Policy>
Mesh<Policy> Mesh<Policy>::cloneWithDofs(const size_type dofPerVertex,
                                         const size_type dofPerElement) const {
  auto mesh = Mesh(_totalNumberOfElements, _totalNumberOfVertices);
  if (_mesh) {
    mesh._mesh = [](const DM &dm) {
      auto clonedDM = DM{};
      Policy::handleError(DMClone(dm, &clonedDM));
      return makeUniqueEntity<Policy>(clonedDM);
    }(_mesh.get());
  }

  auto migration = PetscSF();
  Policy::handleError(DMPlexGetMigrationSF(_mesh.get(), &migration));
  if (migration) {
    auto cloned = PetscSF();
    Policy::handleError(
        PetscSFDuplicate(migration, PETSCSF_DUPLICATE_GRAPH, &cloned));
    Policy::handleError(DMPlexSetMigrationSF(mesh._mesh.get(), cloned));
  }

  mesh.addSection(dofPerVertex, dofPerElement);

  return mesh;
}

template <class Policy>
template <class Container>
void Mesh<Policy>::addChart(const Container &elementVertexIDs,
                            const size_type numberOfVertices, const DM dm) {
  const auto numberOfElements = static_cast<size_type>(elementVertexIDs.size());

  if (Policy::isPrimaryRank()) {
    // create chart
    //
    // element IDs are mapped to mesh IDs:
    // 0, ..., numberOfElements - 1
    //
    // vertex IDs are mapped to mesh IDs:
    // numberOfElements, ..., numberOfElements + numberOfVertices - 1
    Policy::handleError(
        DMPlexSetChart(dm, 0, numberOfVertices + numberOfElements));

    // create data structures
    for (size_type i = 0; i < numberOfElements; ++i) {
      const auto &IDs = elementVertexIDs[i];
      auto numberOfElementVertices = static_cast<size_type>(IDs.size());
      for (const auto id : IDs) {
        if (id == static_cast<decltype(id)>(IGNORE_VERTEX_INDEX))
          numberOfElementVertices--;
      }
      Policy::handleError(DMPlexSetConeSize(dm, i, numberOfElementVertices));
    }
    for (size_type i = numberOfElements;
         i < numberOfElements + numberOfVertices; ++i) {
      Policy::handleError(DMPlexSetConeSize(dm, i, 0));
    }
  }

  Policy::handleError(DMSetUp(dm));

  if (Policy::isPrimaryRank()) {
    // copy connectivity to PETSc mesh
    for (size_type i = 0; i < numberOfElements; ++i) {
      const auto &IDs = elementVertexIDs[i];
      const auto vertexMeshIDs = [&IDs, numberOfElements]() {
        std::vector<size_type> vertexMeshIDs;
        vertexMeshIDs.reserve(IDs.size());
        for (const auto id : IDs) {
          if (id == static_cast<decltype(id)>(IGNORE_VERTEX_INDEX))
            continue;

          assert(id >= 0 &&
                 "Valid indices must be greater than or equal to zero.");
          vertexMeshIDs.push_back(id + numberOfElements);
        }
        return vertexMeshIDs;
      }();
      Policy::handleError(DMPlexSetCone(dm, i, vertexMeshIDs.data()));
    }

    Policy::handleError(DMPlexSymmetrize(dm));
  }
}

template <class Policy>
void Mesh<Policy>::addSection(const size_type dofPerVertex,
                              const size_type dofPerElement) {
  const size_type dimension = topologicalDimension();
  assert(dimension > 0);

  auto section = PetscSection();
  const std::array<size_type, 1> numberOfComponents = {{1}};
  std::vector<size_type> numberOfDofsPerDim(dimension + 1);
  numberOfDofsPerDim.front() = dofPerVertex;
  numberOfDofsPerDim.back() = dofPerElement;

  Policy::handleError(DMPlexCreateSection(_mesh.get(), dimension, 1,
                                          numberOfComponents.data(),
                                          numberOfDofsPerDim.data(), 0, nullptr,
                                          nullptr, nullptr, nullptr, &section));

  Policy::handleError(DMSetDefaultSection(_mesh.get(), section));

  // create default global section
  Policy::handleError(DMGetDefaultGlobalSection(_mesh.get(), &section));
}

template <class Policy> void Mesh<Policy>::distributeMesh(Mesh *const mesh) {
  auto dm = DM();
  auto sf = PetscSF();

  setAdjacencyRules(mesh->_mesh.get());
  Policy::handleError(DMPlexDistribute(mesh->_mesh.get(), 0, &sf, &dm));
  if (dm) {
    Policy::handleError(DMPlexSetMigrationSF(dm, sf));
    mesh->_mesh.reset(dm);
  }
}

template <class Policy> void Mesh<Policy>::setAdjacencyRules(const DM dm) {
  Policy::handleError(DMPlexSetAdjacencyUseCone(dm, PETSC_FALSE));
  Policy::handleError(DMPlexSetAdjacencyUseClosure(dm, PETSC_TRUE));
}

template <class Policy>
Mesh<Policy>::Mesh(const size_type totalNumberOfElements,
                   const size_type totalNumberOfVertices)
    : _mesh(createDM()), _totalNumberOfElements(totalNumberOfElements),
      _totalNumberOfVertices(totalNumberOfVertices) {}

template <class Policy> UniqueEntity<DM> Mesh<Policy>::createDM() {
  auto dm = DM();
  Policy::handleError(DMPlexCreate(Policy::communicator(), &dm));
  return makeUniqueEntity<Policy>(dm);
}

template <class Policy>
typename Mesh<Policy>::size_type Mesh<Policy>::totalNumberOfVertices() const {
  return _totalNumberOfVertices;
}

template <class Policy>
typename Mesh<Policy>::size_type Mesh<Policy>::totalNumberOfElements() const {
  return _totalNumberOfElements;
}

template <class Policy>
typename Mesh<Policy>::size_type Mesh<Policy>::localNumberOfElements() const {
  auto start = size_type{0};
  auto end = size_type{0};
  Policy::handleError(DMPlexGetHeightStratum(_mesh.get(), 0, &start, &end));
  return end - start;
}

template <class Policy>
typename Mesh<Policy>::TopologicalDimension
Mesh<Policy>::topologicalDimension() const {
  auto dim = size_type{};
  Policy::handleError(DMGetDimension(_mesh.get(), &dim));
  return TopologicalDimension(dim);
}

template <class Policy>
typename Mesh<Policy>::CoordinateDimension
Mesh<Policy>::coordinateDimension() const {
  auto dim = size_type{};
  Policy::handleError(DMGetCoordinateDim(_mesh.get(), &dim));
  return CoordinateDimension(dim);
}

template <class Policy>
typename Mesh<Policy>::size_type Mesh<Policy>::localNumberOfVertices() const {
  auto start = size_type{0};
  auto end = size_type{0};
  Policy::handleError(DMPlexGetDepthStratum(_mesh.get(), 0, &start, &end));
  return end - start;
}

template <class Policy> void Mesh<Policy>::print() const {
  Policy::handleError(DMView(_mesh.get(), nullptr));
}

template <class Policy> DM Mesh<Policy>::data() const { return _mesh.get(); }

template <class Policy>
void Mesh<Policy>::copyToLocalVector(
    const distributed<vector_type> &globalVector,
    local<vector_type> *const localVector) const {
  assertCorrectBaseMesh(globalVector.unwrap());
  assertCorrectBaseMesh(localVector->unwrap());

  Policy::handleError(
      DMGlobalToLocalBegin(_mesh.get(), globalVector.unwrap().data(),
                           INSERT_VALUES, localVector->unwrap().data()));
  Policy::handleError(
      DMGlobalToLocalEnd(_mesh.get(), globalVector.unwrap().data(),
                         INSERT_VALUES, localVector->unwrap().data()));
}

template <class Policy>
void Mesh<Policy>::addToGlobalVector(
    const local<vector_type> &localVector,
    distributed<vector_type> *const globalVector) const {
  assertCorrectBaseMesh(localVector.unwrap());
  assertCorrectBaseMesh(globalVector->unwrap());

  Policy::handleError(
      DMLocalToGlobalBegin(_mesh.get(), localVector.unwrap().data(), ADD_VALUES,
                           globalVector->unwrap().data()));
  Policy::handleError(
      DMLocalToGlobalEnd(_mesh.get(), localVector.unwrap().data(), ADD_VALUES,
                         globalVector->unwrap().data()));
}

template <class Policy>
void Mesh<Policy>::copyToGlobalVector(
    const local<vector_type> &localVector,
    distributed<vector_type> *const globalVector) const {
  assertCorrectBaseMesh(localVector.unwrap());
  assertCorrectBaseMesh(globalVector->unwrap());

  Policy::handleError(
      DMLocalToGlobalBegin(_mesh.get(), localVector.unwrap().data(),
                           INSERT_VALUES, globalVector->unwrap().data()));
  Policy::handleError(
      DMLocalToGlobalEnd(_mesh.get(), localVector.unwrap().data(),
                         INSERT_VALUES, globalVector->unwrap().data()));
}

template <class Policy>
typename Mesh<Policy>::size_type
Mesh<Policy>::elementPointIndexToGlobalIndex(const size_type pointIndex) const {
  auto globalIndex = size_type{0};
  auto migration = PetscSF();
  Policy::handleError(DMPlexGetMigrationSF(_mesh.get(), &migration));
  if (migration) {
    const PetscSFNode *nodes = nullptr;
    Policy::handleError(
        PetscSFGetGraph(migration, nullptr, nullptr, nullptr, &nodes));
    return nodes[pointIndex].index;
  }

  return pointIndex;
}

template <class Policy>
typename Mesh<Policy>::size_type
Mesh<Policy>::vertexPointIndexToGlobalIndex(const size_type pointIndex) const {
  return elementPointIndexToGlobalIndex(pointIndex) - _totalNumberOfElements;
}

template <class Policy>
typename Mesh<Policy>::size_type
Mesh<Policy>::numberOfVertices(const size_type elementPointIndex) const {
  size_type size = 0;
  Policy::handleError(DMPlexGetConeSize(_mesh.get(), elementPointIndex, &size));
  return size;
}

template <class Policy>
std::vector<typename Mesh<Policy>::size_type>
Mesh<Policy>::vertices(const size_type elementPointIndex) const {
  const auto size = numberOfVertices(elementPointIndex);
  const size_type *cone = nullptr;
  Policy::handleError(DMPlexGetCone(_mesh.get(), elementPointIndex, &cone));
  std::vector<size_type> output(size);
  std::transform(cone, cone + size, output.begin(),
                 [this](const size_type entityPointIndex) {
                   return vertexPointIndexToGlobalIndex(entityPointIndex);
                 });
  return output;
}

template <class Policy>
std::pair<typename Mesh<Policy>::size_type, typename Mesh<Policy>::size_type>
Mesh<Policy>::localDofLineRange(const size_type entityPointIndex) const {
  auto output = std::pair<size_type, size_type>();
  Policy::handleError(DMPlexGetPointLocal(_mesh.get(), entityPointIndex,
                                          &output.first, &output.second));
  return output;
}

template <class Policy>
std::pair<typename Mesh<Policy>::size_type, typename Mesh<Policy>::size_type>
Mesh<Policy>::globalDofLineRange(const size_type entityPointIndex) const {
  auto output = std::pair<size_type, size_type>();
  Policy::handleError(DMPlexGetPointGlobal(_mesh.get(), entityPointIndex,
                                           &output.first, &output.second));
  return output;
}

template <class Policy>
typename Mesh<Policy>::size_type
Mesh<Policy>::numberOfDofs(const size_type entityPointIndex) const {
  auto output = size_type{0};
  auto section = PetscSection();
  Policy::handleError(DMGetDefaultSection(_mesh.get(), &section));

  Policy::handleError(PetscSectionGetDof(section, entityPointIndex, &output));

  return output;
}

template <class Policy>
void Mesh<Policy>::copyEntityData(const size_type entityPointIndex,
                                  const local<vector_type> &localVector,
                                  std::vector<value_type> *const data) const {
  assertCorrectBaseMesh(localVector.unwrap());
  auto size = size_type{0};
  Policy::handleError(DMPlexVecGetClosure(_mesh.get(), nullptr,
                                          localVector.unwrap().data(),
                                          entityPointIndex, &size, nullptr));
  data->resize(size);
  if (size == 0)
    return;

  auto dataPtr = data->data();
  Policy::handleError(DMPlexVecGetClosure(_mesh.get(), nullptr,
                                          localVector.unwrap().data(),
                                          entityPointIndex, &size, &dataPtr));
}

template <class Policy>
void Mesh<Policy>::addEntityData(const size_type entityPointIndex,
                                 const std::vector<value_type> &data,
                                 local<vector_type> *const localVector) const {
  assertCorrectBaseMesh(localVector->unwrap());
  Policy::handleError(
      DMPlexVecSetClosure(_mesh.get(), nullptr, localVector->unwrap().data(),
                          entityPointIndex, data.data(), ADD_VALUES));
}

template <class Policy>
void Mesh<Policy>::setEntityData(const size_type entityPointIndex,
                                 const std::vector<value_type> &data,
                                 local<vector_type> *const localVector) const {
  assertCorrectBaseMesh(localVector->unwrap());
  Policy::handleError(
      DMPlexVecSetClosure(_mesh.get(), nullptr, localVector->unwrap().data(),
                          entityPointIndex, data.data(), INSERT_VALUES));
}

template <class Policy>
void Mesh<Policy>::addEntityMatrix(const size_type entityPointIndex,
                                   const std::vector<value_type> &data,
                                   matrix_type *const matrix) const {
  Policy::handleError(DMPlexMatSetClosure(_mesh.get(), nullptr, nullptr,
                                          matrix->data(), entityPointIndex,
                                          data.data(), ADD_VALUES));
}

template <class Policy>
void Mesh<Policy>::assertCorrectBaseMesh(const vector_type &vector) const {
#ifndef NDEBUG
  const auto extractDM = [](const vector_type &vector) {
    auto dm = DM{};
    Policy::handleError(VecGetDM(vector.data(), &dm));
    return dm;
  };
  const auto dm = extractDM(vector);
  assert((!dm || dm == _mesh.get()) &&
         "The vector was created using a different mesh.");
#else
  static_cast<void>(vector);
#endif
}
} // namespace cpppetsc
} // namespace ae108