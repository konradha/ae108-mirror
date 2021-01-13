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

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/elements/mesh/Mesh.h"
#include "ae108/elements/mesh/generate_triangle_mesh.h"
#include "ae108/elements/tensor/Tensor.h"
#include <array>

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Mesh = cpppetsc::Mesh<Policy>;

// In this example we will mesh a rectangle with triangles and distribute the
// result with PETSc.

// Let's specify the parameters of this mesh.

constexpr auto dimension = Mesh::size_type{2};
constexpr auto dof_per_vertex = Mesh::size_type{2};

// We choose the rectangle to have x-length=3 and y-length=2.

constexpr auto rectangle_size = elements::tensor::Tensor<double, 2>{{3., 2.}};

// Also, we want to generate 2 * 4 triangles in x-direction, and 2 * 1 triangles
// in y-direction.

constexpr auto mesh_granularity =
    elements::tensor::Tensor<std::size_t, 2>{{4, 1}};

int main(int argc, char **argv) {
  // PETSc must be initialized before using it.
  Policy::handleError(PetscInitialize(&argc, &argv, NULL, NULL));

  // We use a scope around our computation to make sure everything is cleaned up
  // before we call PetscFinalize.
  {
    // First we generate the mesh locally. This mesh will look the same on every
    // rank.

    const auto generated_mesh = elements::mesh::generate_triangle_mesh(
        rectangle_size, mesh_granularity);

    // The number of positions in the mesh provides the number of vertices here.

    const auto number_of_vertices = generated_mesh.number_of_positions();

    // Let's create the distributed PETSc-mesh using the generated connectivity
    // and print it to the console. We expect 2 * 4 2-cells (triangles) and
    // 2 * (4 + 1) 0-cells (vertices). In parallel, some vertices are
    // shared between ranks and we therefore expect a sum of at least 10
    // 0-cells.

    const auto mesh =
        Mesh::fromConnectivity(dimension, generated_mesh.connectivity(),
                               number_of_vertices, dof_per_vertex);
    mesh.print();

    // Before we continue, we wait until all ranks are done printing.

    Policy::handleError(PetscBarrier(nullptr));

    // The PETSc mesh does not contain vertex positions, only the generated mesh
    // does. Therefore, use the generated mesh to access vertex positions using
    // an index. To test that this works, we print the position of the second
    // vertex (index 1) on every rank that shares it. We expect this vertex to
    // be at (0., 2.).

    for (const auto &vertex : mesh.localVertices()) {
      if (vertex.index() != 1) {
        continue;
      }

      // We have found the vertex with index 1 locally.

      const auto position = generated_mesh.position_of_vertex(1);
      printf("The second vertex is at (%f, %f).\n", position.at(0),
             position.at(1));
    }

    // If you choose enough processes (in my case 4 works), you will see that
    // the vertex is shared by two ranks because more than one rank prints the
    // message.
  }

  Policy::handleError(PetscFinalize());
}
