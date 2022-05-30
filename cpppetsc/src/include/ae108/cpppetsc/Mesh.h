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

#include "ae108/cpppetsc/LocalElementView.h"
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
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <utility>
#include <vector>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class Mesh {
public:
  using policy_type = Policy;

  using size_type = PetscInt;
  using value_type = PetscScalar;
  using real_type = PetscReal;
  using vector_type = Vector<Policy>;
  using matrix_type = Matrix<Policy>;

  static constexpr size_type IGNORE_VERTEX_INDEX = -1;

  /**
   * @brief Creates a Mesh from the connectivity representation.
   *
   * @tparam Container Requirements: size(); Container::value_type requirements:
   * size(), data().
   *
   * @param dimension The coordinate dimension.
   *
   * @param elementVertexIDs A container of containers: the container at index i
   * contains the vertex IDs of the element with index i. Each vertex ID must
   * lie in [0, numberOfVertices) u {IGNORE_VERTEX_INDEX}. IDs of value
   * IGNORE_VERTEX_INDEX are ignored.
   *
   * @param numberOfVertices The maximum vertex ID in the mesh.
   */
  template <class Container>
  static Mesh fromConnectivity(const size_type dimension,
                               const Container &elementVertexIDs,
                               const size_type numberOfVertices,
                               const size_type dofPerVertex,
                               const size_type dofPerElement = 0);

  /**
   * @brief Clone the mesh with a different default section.
   */
  Mesh cloneWithDofs(const size_type dofPerVertex,
                     const size_type dofPerElement) const;

  /**
   * @brief Makes it possible to iterate over (views of) the local elements.
   */
  auto localElements() const;

  /**
   * @brief Makes it possible to iterate over (views of) the local vertices.
   */
  auto localVertices() const;

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
   * @brief Returns the dimension of the coordinates.
   */
  size_type coordinateDimension() const;

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
   * @brief Returns a global vector in canonical row order (the degrees of
   * freedom of the elements sorted by global index, followed by the degrees of
   * freedom of the vertices sorted by global index) transformed from the
   * provided vector in PETSc's row order.
   */
  distributed<vector_type>
  toCanonicalOrder(const distributed<vector_type> &vector) const;

  /**
   * @brief Returns a matrix in canonical row/column order (the degrees of
   * freedom of the elements sorted by global index, followed by the degrees of
   * freedom of the vertices sorted by global index) transformed from the
   * provided matrix in PETSc's row/column order.
   */
  matrix_type toCanonicalOrder(const matrix_type &matrix) const;

  /**
   * @brief Returns a global vector in PETSc's row order transformed from the
   * provided vector in canonical row order (the degrees of freedom of the
   * elements sorted by global index, followed by the degrees of freedom of the
   * vertices sorted by global index).
   */
  distributed<vector_type>
  fromCanonicalOrder(const distributed<vector_type> &vector) const;

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
   * @brief Distributes the provided mesh over the MPI nodes.
   */
  static void distributeMesh(Mesh *const);

  /**
   * @brief Returns the layout of a global vector associated with this mesh.
   */
  PetscLayout globalVectorLayout() const;

  /**
   * @brief Returns a vector that contains the canonical row index for every
   * locally-owned row in a global vector.
   */
  std::vector<size_type>
  canonicalRowIndices(const size_type dofPerVertex,
                      const size_type dofPerElement) const;

  /**
   * @brief Computes the transform between PETSc's row order and the order
   * specified by `targetRows` for locally-owned rows.
   */
  UniqueEntity<PetscSF>
  createReorderingSF(const std::vector<size_type> &targetRows) const;

  /**
   * @brief Stores `globalToNatural` as the "GlobalToNaturalSF" of the
   * DM.
   */
  void setGlobalToNaturalSF(UniqueEntity<PetscSF> globalToNatural);

  /**
   * @brief Asserts that the vector has been created with this mesh.
   */
  void assertCorrectBaseMesh(const vector_type &vector) const;

  UniqueEntity<DM> _mesh;
  UniqueEntity<PetscSF> _migration;
  UniqueEntity<PetscSection> _section;
  UniqueEntity<PetscSF> _globalToNatural;
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
#include "ae108/cpppetsc/MeshDataProvider.h"
#include "ae108/cpppetsc/Vector.h"

namespace ae108 {
namespace cpppetsc {

template <class Policy> auto Mesh<Policy>::localElements() const {
  auto start = size_type{0};
  auto stop = size_type{0};
  Policy::handleError(DMPlexGetHeightStratum(_mesh.get(), 0, &start, &stop));

  namespace rv = ranges::cpp20::views;
  const auto provider = createDataProviderFromMesh(this);
  return rv::iota(start, stop) | rv::transform([provider](const size_type id) {
           return LocalElementView<Mesh>(provider, id);
         });
}

template <class Policy> auto Mesh<Policy>::localVertices() const {
  auto start = size_type{0};
  auto stop = size_type{0};
  Policy::handleError(DMPlexGetDepthStratum(_mesh.get(), 0, &start, &stop));

  namespace rv = ranges::cpp20::views;
  const auto provider = createDataProviderFromMesh(this);
  return rv::iota(start, stop) | rv::transform([provider](const size_type id) {
           return LocalVertexView<Mesh>(provider, id);
         });
}

template <class Policy>
template <class Container>
Mesh<Policy> Mesh<Policy>::fromConnectivity(const size_type dimension,
                                            const Container &elementVertexIDs,
                                            const size_type numberOfVertices,
                                            const size_type dofPerVertex,
                                            const size_type dofPerElement) {
  auto mesh =
      Mesh(static_cast<size_type>(elementVertexIDs.size()), numberOfVertices);

  Policy::handleError(DMSetCoordinateDim(mesh._mesh.get(), dimension));
  // set dimension to 1 to use two levels in the graph:
  // the vertices and the elements
  Policy::handleError(DMSetDimension(mesh._mesh.get(), 1));

  addChart(elementVertexIDs, numberOfVertices, mesh._mesh.get());

  // calculate the strata
  Policy::handleError(DMPlexStratify(mesh._mesh.get()));

  distributeMesh(&mesh);

  mesh.addSection(dofPerVertex, dofPerElement);

  mesh.setGlobalToNaturalSF(mesh.createReorderingSF(
      mesh.canonicalRowIndices(dofPerVertex, dofPerElement)));

  return mesh;
}

template <class Policy> PetscLayout Mesh<Policy>::globalVectorLayout() const {
  const auto vector = [this] {
    auto vec = Vec();
    Policy::handleError(DMGetGlobalVector(_mesh.get(), &vec));
    return UniqueEntity<Vec>(vec, [&](Vec vec) {
      Policy::handleError(DMRestoreGlobalVector(_mesh.get(), &vec));
    });
  }();

  auto layout = PetscLayout();
  Policy::handleError(VecGetLayout(vector.get(), &layout));
  return layout;
}

template <class Policy>
std::vector<typename Mesh<Policy>::size_type>
Mesh<Policy>::canonicalRowIndices(const size_type dofPerVertex,
                                  const size_type dofPerElement) const {
  using Range = std::pair<size_type, size_type>;

  const auto localRange = [this] {
    auto result = Range();
    Policy::handleError(PetscLayoutGetRange(globalVectorLayout(), &result.first,
                                            &result.second));
    return result;
  }();

  auto result = std::vector<size_type>(localRange.second - localRange.first,
                                       size_type{-1});

  const auto insert = [&](const Range &range, const size_type offset) {
    for (auto row = std::max(range.first, localRange.first);
         row < std::min(range.second, localRange.second); ++row) {
      const auto dof = row - range.first;
      result[row - localRange.first] = dof + offset;
    }
  };

  for (const auto &element : localElements()) {
    insert(element.globalDofLineRange(), element.index() * dofPerElement);
  }
  for (const auto &vertex : localVertices()) {
    insert(vertex.globalDofLineRange(), _totalNumberOfElements * dofPerElement +
                                            vertex.index() * dofPerVertex);
  }

  return result;
}

template <class Policy>
UniqueEntity<PetscSF> Mesh<Policy>::createReorderingSF(
    const std::vector<size_type> &targetRows) const {
  const auto inverse = []() {
    auto result = PetscSF();
    Policy::handleError(PetscSFCreate(Policy::communicator(), &result));
    return makeUniqueEntity<Policy>(result);
  }();
  Policy::handleError(PetscSFSetGraphLayout(
      inverse.get(), globalVectorLayout(), targetRows.size(), nullptr,
      PETSC_USE_POINTER, targetRows.data()));
  Policy::handleError(PetscSFSetUp(inverse.get()));

  auto result = PetscSF();
  Policy::handleError(PetscSFCreateInverseSF(inverse.get(), &result));
  return makeUniqueEntity<Policy>(result);
}

template <class Policy>
void Mesh<Policy>::setGlobalToNaturalSF(UniqueEntity<PetscSF> globalToNatural) {
  _globalToNatural = std::move(globalToNatural);
  Policy::handleError(
      DMPlexSetGlobalToNaturalSF(_mesh.get(), _globalToNatural.get()));
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
    mesh._migration = [&]() {
      auto sf = PetscSF();
      Policy::handleError(
          PetscSFDuplicate(migration, PETSCSF_DUPLICATE_GRAPH, &sf));
      return makeUniqueEntity<Policy>(sf);
    }();
    Policy::handleError(
        DMPlexSetMigrationSF(mesh._mesh.get(), mesh._migration.get()));
  }

  mesh.addSection(dofPerVertex, dofPerElement);

  mesh.setGlobalToNaturalSF(mesh.createReorderingSF(
      mesh.canonicalRowIndices(dofPerVertex, dofPerElement)));

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
  const std::array<size_type, 1> numberOfComponents = {{1}};
  std::vector<size_type> numberOfDofsPerDim(2);
  numberOfDofsPerDim.front() = dofPerVertex;
  numberOfDofsPerDim.back() = dofPerElement;

  Policy::handleError(DMSetNumFields(_mesh.get(), numberOfComponents.size()));

  _section = [&]() {
    auto s = PetscSection();
    Policy::handleError(DMPlexCreateSection(
        _mesh.get(), nullptr /* label */, numberOfComponents.data(),
        numberOfDofsPerDim.data(), 0 /* number of boundary conditions */,
        nullptr /* boundary conditions */,
        nullptr /* boundary condition components */,
        nullptr /* boundary condition points */, nullptr /* permutation */,
        &s));
    return makeUniqueEntity<Policy>(s);
  }();

  Policy::handleError(DMSetSection(_mesh.get(), _section.get()));
}

template <class Policy> void Mesh<Policy>::distributeMesh(Mesh *const mesh) {
  auto dm = DM();
  auto sf = PetscSF();

  Policy::handleError(DMPlexDistribute(mesh->_mesh.get(), 0, &sf, &dm));
  if (dm) {
    Policy::handleError(DMPlexSetMigrationSF(dm, sf));
    mesh->_mesh.reset(dm);
  }
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
typename Mesh<Policy>::size_type Mesh<Policy>::coordinateDimension() const {
  auto dim = size_type{};
  Policy::handleError(DMGetCoordinateDim(_mesh.get(), &dim));
  return dim;
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
distributed<typename Mesh<Policy>::vector_type>
Mesh<Policy>::toCanonicalOrder(const distributed<vector_type> &vector) const {
  auto globalToNatural = PetscSF{};
  Policy::handleError(
      DMPlexGetGlobalToNaturalSF(_mesh.get(), &globalToNatural));

  auto result = tag<DistributedTag>(Vector<Policy>::fromLayoutOf(vector));
  Policy::handleError(DMPlexGlobalToNaturalBegin(
      _mesh.get(), vector.unwrap().data(), result.unwrap().data()));
  Policy::handleError(DMPlexGlobalToNaturalEnd(
      _mesh.get(), vector.unwrap().data(), result.unwrap().data()));
  return result;
}

template <class Policy>
typename Mesh<Policy>::matrix_type
Mesh<Policy>::toCanonicalOrder(const matrix_type &matrix) const {
  const auto mapping = [&]() {
    auto sf = PetscSF{};
    Policy::handleError(DMPlexGetGlobalToNaturalSF(_mesh.get(), &sf));
    auto mapping = ISLocalToGlobalMapping();
    Policy::handleError(ISLocalToGlobalMappingCreateSF(
        sf, matrix.localRowRange().first, &mapping));
    return makeUniqueEntity<Policy>(mapping);
  }();

  const auto reorder = [&](const auto &in) {
    auto out = IS();
    Policy::handleError(
        ISLocalToGlobalMappingApplyIS(mapping.get(), in.get(), &out));
    return makeUniqueEntity<Policy>(out);
  };

  const auto indices = [&](const auto n) {
    auto is = IS();
    Policy::handleError(ISCreateStride(Policy::communicator(), n, 0, 1, &is));
    return reorder(makeUniqueEntity<Policy>(is));
  };

  const auto rows = indices(matrix.localSize().first);
  const auto cols = indices(matrix.localSize().second);

  auto mat = Mat();
  Policy::handleError(MatCreateSubMatrix(matrix.data(), rows.get(), cols.get(),
                                         MAT_INITIAL_MATRIX, &mat));
  return matrix_type(makeUniqueEntity<Policy>(mat));
}

template <class Policy>
distributed<typename Mesh<Policy>::vector_type>
Mesh<Policy>::fromCanonicalOrder(const distributed<vector_type> &vector) const {
  auto globalToNatural = PetscSF{};
  Policy::handleError(
      DMPlexGetGlobalToNaturalSF(_mesh.get(), &globalToNatural));

  auto result = tag<DistributedTag>(Vector<Policy>::fromLayoutOf(vector));
  Policy::handleError(DMPlexNaturalToGlobalBegin(
      _mesh.get(), vector.unwrap().data(), result.unwrap().data()));
  Policy::handleError(DMPlexNaturalToGlobalEnd(
      _mesh.get(), vector.unwrap().data(), result.unwrap().data()));
  return result;
}

template <class Policy>
void Mesh<Policy>::assertCorrectBaseMesh(const vector_type &vector) const {
#if !defined(NDEBUG) && !defined(AE108_CPPPETSC_DISABLE_BASE_MESH_CHECK)
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
