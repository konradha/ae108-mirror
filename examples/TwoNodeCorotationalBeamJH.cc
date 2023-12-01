// © 2023 ETH Zurich, Mechanics and Materials Lab
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

// TODO: FINISH
// currently just copying most of TimoshenkoBeam as an example

#include "ae108/assembly/Assembler.h"
#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/Viewer.h"
#include "ae108/cpppetsc/createVectorFromSource.h"
#include "ae108/cpppetsc/setName.h"
#include "ae108/cpppetsc/writeToViewer.h"
#include "ae108/elements/TwoNodeCorotationalBeamJHElement.h"
#include "ae108/elements/tensor/as_vector.h"
#include "ae108/solve/NonlinearSolver.h"


using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using Viewer = cpppetsc::Viewer<Policy>;
using BoundaryCondition = cpppetsc::MeshBoundaryCondition<Mesh>;

// In this example we will simulate six elements, A, B, C, D, E, F.
// Each of these elements has two vertices.

//  4 --D---3
//    \     |  |
//      F   E    C
//        \ |      |
//  0---A---1---B---2

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{2};
constexpr auto number_of_elements = Mesh::size_type{6};
constexpr auto number_of_vertices = Mesh::size_type{5};
constexpr auto coordinate_dimension = Mesh::size_type{2};
constexpr auto dof_per_vertex = Mesh::size_type{3};
constexpr auto dof_per_element = Mesh::size_type{0};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1}}, // vertices of element A
    {{1, 2}}, // vertices of element B
    {{2, 3}}, // vertices of element C
    {{3, 4}}, // vertices of element D
    {{1, 3}}, // vertices of element E
    {{1, 4}}, // vertices of element F
}};

// Each of the vertices is located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::real_type, coordinate_dimension>,
               number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0., 0.}},
    {{1., 0.}},
    {{2., 0.}},
    {{1., 1.}},
    {{0., 1.}},
}};

using Element =
    elements::TwoNodeCorotationalBeamJHElement<coordinate_dimension, Vector::value_type,
                                    Vector::real_type>;

using Properties =
    elements::TwoNodeCorotationalBeamJHProperties<Mesh::real_type, coordinate_dimension>;

constexpr Mesh::real_type young_modulus = 1.;
constexpr Mesh::real_type poisson_ratio = 0.3;
constexpr Mesh::real_type shear_modulus =
    young_modulus / (2 * (1 + poisson_ratio));
constexpr Mesh::real_type shear_correction_factor_y = 1.2;
constexpr Mesh::real_type thickness = 1.;
constexpr Mesh::real_type width = 0.1;
constexpr Mesh::real_type area = width * thickness;
constexpr Mesh::real_type area_moment_z =
    thickness * width * width * width / 12.;

using Assembler =
    assembly::Assembler<Element, assembly::DefaultFeaturePlugins, Policy>;


using Solver = solve::NonlinearSolver<Assembler>;

int main() {
  return 0;
}