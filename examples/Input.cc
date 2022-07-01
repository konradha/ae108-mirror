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
#include "ae108/cpppetsc/readFromViewer.h"
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

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // First we create a mesh.
  const auto mesh = Mesh::fromConnectivity(coordinate_dimension, connectivity,
                                           number_of_vertices, dof_per_vertex,
                                           dof_per_element);

  // Let's create a global vector and fill it with the vertex indices.
  using DataSource = std::function<void(Mesh::size_type, Mesh::value_type *)>;
  auto data = cpppetsc::createVectorFromSource(
      mesh, dof_per_vertex,
      DataSource([&](const Mesh::size_type index, Mesh::value_type *const out) {
        std::fill_n(out, dof_per_vertex, static_cast<Mesh::value_type>(index));
      }));
  cpppetsc::setName("vertex_index_data", &data);

  // To be able to compare it with the read vector we print it to the console.
  data.unwrap().print();

  // We write this vector to "input.ae108".
  {
    auto viewer = Viewer::fromHdf5FilePath("input.ae108", Viewer::Mode::write);
    cpppetsc::writeToViewer(data, &viewer);
  }

  // Finally we read this vector from the file to a new vector called `input`.
  auto input = Vector::fromGlobalMesh(mesh);
  {
    const auto viewer =
        Viewer::fromHdf5FilePath("input.ae108", Viewer::Mode::read);
    cpppetsc::setName("vertex_index_data", &input);
    cpppetsc::readFromViewer(viewer, &input);
  }

  // We expect the same vector as printed above when printing the result.
  input.unwrap().print();
}