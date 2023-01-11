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

#include "ae108/meshing/cppgmsh/extract_mesh.h"
#include "ae108/meshing/cppgmsh/get_elements_in.h"
#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include <gmsh.h>
#include <range/v3/view/enumerate.hpp>
namespace ae108 {
namespace meshing {
namespace cppgmsh {

template <std::size_t coordinate_dimension>
std::tuple<std::vector<std::array<double, coordinate_dimension>>,
           std::vector<std::vector<std::size_t>>,
           std::map<std::size_t, std::size_t>,
           std::map<std::size_t, std::size_t>>
extract_mesh(const std::size_t element_dimension) noexcept {

  const auto [nodeTagToIndex,
              positions] = [](const std::size_t element_dimension) {
    auto nodes = get_nodes_of<coordinate_dimension>({element_dimension, -1});

    auto nodeTagToIndex = std::map<std::size_t, std::size_t>();
    auto positions =
        std::vector<std::array<double, coordinate_dimension>>(nodes.size());
    for (const auto &&[index, node] : ranges::views::enumerate(nodes)) {
      positions[index] = node.position;
      nodeTagToIndex[node.id] = index;
    }

    return std::make_tuple(nodeTagToIndex, positions);
  }(element_dimension);

  const auto [elementTagToIndex, connectivity] =
      [&nodeTagToIndex](const std::size_t element_dimension) {
        const auto elements = get_elements_in({element_dimension, -1});

        auto elementTagToIndex = std::map<std::size_t, std::size_t>();
        auto connectivity = std::vector<std::vector<std::size_t>>();
        connectivity.reserve(elements.size());
        for (auto element : elements) {
          std::vector<std::size_t> nodeTags;
          gmsh::model::mesh::getElement(element.second, element.first,
                                        nodeTags);

          std::vector<std::size_t> nodeIndices(nodeTags.size());
          for (auto &&[index, nodeTag] : ranges::views::enumerate(nodeTags))
            nodeIndices[index] = nodeTagToIndex.at(nodeTag);

          elementTagToIndex[element.second] = connectivity.size();
          connectivity.push_back(nodeIndices);
        }

        return std::make_tuple(elementTagToIndex, connectivity);
      }(element_dimension);

  return {positions, connectivity, nodeTagToIndex, elementTagToIndex};
}

template std::tuple<
    std::vector<std::array<double, 3>>, std::vector<std::vector<std::size_t>>,
    std::map<std::size_t, std::size_t>, std::map<std::size_t, std::size_t>>
extract_mesh<3>(const std::size_t element_dimension) noexcept;

template std::tuple<
    std::vector<std::array<double, 2>>, std::vector<std::vector<std::size_t>>,
    std::map<std::size_t, std::size_t>, std::map<std::size_t, std::size_t>>
extract_mesh<2>(const std::size_t element_dimension) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108