// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/meshing/bounding_box_of.h"
#include "ae108/meshing/construct_periodic_point_cloud.h"
#include "ae108/meshing/construct_voronoi_cell.h"
#include "ae108/meshing/cppgmsh/Context.h"
#include "ae108/meshing/cppgmsh/construct_cylinders.h"
#include "ae108/meshing/cppgmsh/construct_polyhedron.h"
#include "ae108/meshing/cppgmsh/extract_mesh.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_entities_in.h"
#include "ae108/meshing/cppgmsh/get_periodic_nodes_of.h"
#include "ae108/meshing/cppgmsh/intersect_entities.h"
#include "ae108/meshing/cppgmsh/set_domain_entities_periodic.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include "ae108/meshing/cppgmsh/write_file.h"

#include <iostream>

using namespace ae108::meshing;

// In this example we will generate a periodic mesh of a representative volume
// element as it is often required for homogenization or Bloch wave analysis of
// a periodic medium.

constexpr std::size_t dimension = 3;
using Point = std::array<double, dimension>;

// A 3D periodic medium with translational symmetry is commonly described as a
// Bravais lattice by three translation vectors. Three well-known lattices are
// the simple cubic (sc), the body centered cubic (bcc), and the face centered
// cubic (fcc) lattice. https://en.wikipedia.org/wiki/Bravais_lattice

constexpr auto lattice_vectors =
    std::array<std::array<double, dimension>, dimension>{
        {{0.5, 0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, -0.5}}}; // bcc
//      {{1., 1., 0.}, {0., 1., 1.}, {1., 0., 1.}}}; //fcc
//      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; //sc

int main(int argc, char **argv) {

  // In a first step, we need to initialize the meshing engine "cppgmsh"

  const auto gmshContext = cppgmsh::Context(argc, argv);

  // Let us now define the domain of our representative volume element (rve).
  // For a perodic medium, like above Bravais lattice, the smallest possible
  // unit cell is a primitive unit cell: Two frequent choices are the
  // parallelpiped or the Wigner-Seitz (a type of Voronoi cell) primitive cell.
  // https://en.wikipedia.org/wiki/Unit_cell

  std::vector<std::pair<std::size_t, std::size_t>> periodic_faces;
  const auto domain = construct_voronoi_cell(
      construct_periodic_point_cloud(lattice_vectors, {0, 0, 0}),
      &periodic_faces);
  const auto domain_entity = cppgmsh::construct_polyhedron(domain);

  // Next, we define the solid within the primitive unit cell. We choose to
  // place cylinders on the edges of the domain. The solid may also span several
  // unit cells like in a finite lattice.

  const auto solid_entity = cppgmsh::construct_cylinders(
      domain.vertices, domain.edges,
      std::vector<double>(domain.edges.size(), 0.05), true);

  // Finally, we construct the rve by cropping all solid outside of the
  // domain of the primitive unit cell.

  const auto rve_entity =
      cppgmsh::intersect_entities({solid_entity}, {domain_entity});

  // We sync the constructive solid geometry engine and write a CAD file.

  cppgmsh::synchronize();
  //   cppgmsh::write_file("periodic.step");

  // To obtain a periodic mesh (e.g. for homogenization or Bloch wave analysis),
  // we define the periodic domain faces.

  for (const auto &face_pair : periodic_faces)
    cppgmsh::set_domain_entities_periodic(
        bounding_box_of(domain.vertices_of_face(face_pair.first)),
        bounding_box_of(domain.vertices_of_face(face_pair.second)), 2);

  // Finally, we generate the periodic FE mesh.

  cppgmsh::set_granularity(1. / 10);
  cppgmsh::generate_mesh();

  // We save the mesh directly as "periodic.vtk". Other formats are available.

  cppgmsh::write_file("periodic.vtk");

  // To use the mesh in an ae108 FE application, we extract the positions and
  // connectivities of the mesh elements.

  const auto [positions, connectivity, nodeTagToIndex, elementTagToIndex] =
      cppgmsh::extract_mesh<dimension>();

  for (const auto &face_pair : periodic_faces) {
    std::pair<int, int> numper_of_nodes;
    for (const auto &target_surface : cppgmsh::get_entities_in(
             bounding_box_of(domain.vertices_of_face(face_pair.second)), 2)) {
      const auto nodes = cppgmsh::get_periodic_nodes_of(target_surface);
      assert(nodes.first.size() == nodes.second.size());
      numper_of_nodes.first += nodes.first.size();
      numper_of_nodes.second += nodes.second.size();
    }
    std::cout << "# target face (" << face_pair.second
              << ") nodes: " << numper_of_nodes.second << std::endl
              << "# source face (" << face_pair.first
              << ") nodes: " << numper_of_nodes.first << std::endl
              << std::endl;
  }
}
