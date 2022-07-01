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

#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/Viewer.h"
#include "ae108/cpppetsc/createVectorFromSource.h"
#include "ae108/cpppetsc/setName.h"
#include "ae108/cpppetsc/writeToViewer.h"
#include <algorithm>
#include <array>

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using Viewer = cpppetsc::Viewer<Policy>;

// In this example we will save data for two elements, A and B, to a HDF5 file.
// Each of these elements has four vertices.

//  3 ------- 2 ------- 5
//  |         |         |
//  |    A    |    B    |
//  0---------1---------4

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{4};
constexpr auto number_of_elements = Mesh::size_type{2};
constexpr auto number_of_vertices = Mesh::size_type{6};
constexpr auto coordinate_dimension = Mesh::size_type{3};
constexpr auto dof_per_vertex = Mesh::size_type{1};
constexpr auto dof_per_element = Mesh::size_type{0};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1, 2, 3}}, // vertices of element A
    {{1, 4, 5, 2}}, // vertices of element B
}};

// Each of the vertices is located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::value_type, coordinate_dimension>,
               number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0., 0., 0.}},
    {{1., 0., 0.}},
    {{1., 1., 0.}},
    {{0., 1., 0.}},
    {{2., 0., 0.}},
    {{2., 1., 0.}},
}};

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // First we create a mesh.
  const auto mesh = Mesh::fromConnectivity(coordinate_dimension, connectivity,
                                           number_of_vertices, dof_per_vertex,
                                           dof_per_element);

  // Then we create a vector of coordinates from `vertex_positions`.
  using DataSource = std::function<void(Mesh::size_type, Mesh::value_type *)>;
  auto coordinates = cpppetsc::createVectorFromSource(
      mesh, coordinate_dimension,
      DataSource([&](const Mesh::size_type index, Mesh::value_type *const out) {
        const auto &position = vertex_positions.at(index);
        std::copy(position.begin(), position.end(), out);
      }));
  cpppetsc::setName("coordinates", &coordinates);

  // Let's create a global vector and fill it with the vertex indices.
  constexpr auto data_per_vertex = coordinate_dimension + 2;
  auto data = cpppetsc::createVectorFromSource(
      mesh, data_per_vertex,
      DataSource([&](const Mesh::size_type index, Mesh::value_type *const out) {
        std::fill_n(out, data_per_vertex, static_cast<Mesh::value_type>(index));
      }));
  cpppetsc::setName("vertex_index_data", &data);

  // Finally we write the results to "output.ae108".
  // We start with writing the mesh.
  auto viewer = Viewer::fromHdf5FilePath("output.ae108", Viewer::Mode::write);
  cpppetsc::writeToViewer(mesh, coordinates, &viewer);

  // Then we write two vector fields of different dimension.
  cpppetsc::writeToViewer(data, &viewer);        // dimension 5
  cpppetsc::writeToViewer(coordinates, &viewer); // dimension 3

  // Note that Paraview supports both "point arrays" and "cell arrays"
  auto point_array =
      Mesh::vector_type::fromGlobalMesh(mesh.cloneWithDofs(1, 0));
  cpppetsc::setName("point_array", &point_array);
  cpppetsc::writeToViewer(point_array, &viewer);
  auto cell_array = Mesh::vector_type::fromGlobalMesh(mesh.cloneWithDofs(0, 1));
  cpppetsc::setName("cell_array", &cell_array);
  cpppetsc::writeToViewer(cell_array, &viewer);

  fprintf(stderr, "The data has been written to the file \"output.ae108\".\n");
}
