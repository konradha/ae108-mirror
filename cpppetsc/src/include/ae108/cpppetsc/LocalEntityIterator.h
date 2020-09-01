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

#include <iterator>
#include <vector>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief A bidirectional iterator to traverse the local entities of a mesh.
 * Dereferencing the iterator yields a view of a local entity.
 *
 * @remark Changes to the mesh invalidate the iterator.
 */
template <class MeshType, template <class> class ViewTemplate>
class LocalEntityIterator {
public:
  using mesh_type = MeshType;

  using View = ViewTemplate<LocalEntityIterator>;

  using value_type = View;
  using iterator_category = std::bidirectional_iterator_tag;
  using difference_type = typename mesh_type::size_type;
  using reference = const View &;
  using pointer = const View *;

  LocalEntityIterator() = default;
  explicit LocalEntityIterator(
      const MeshType *const mesh,
      const typename mesh_type::size_type elementPointIndex);

  LocalEntityIterator &operator++();

  LocalEntityIterator operator++(int);

  LocalEntityIterator &operator--();

  LocalEntityIterator operator--(int);

  reference operator*() const;

  pointer operator->() const;

  friend inline bool operator==(const LocalEntityIterator &lhs,
                                const LocalEntityIterator &rhs) {
    return lhs._view == rhs._view;
  }

  friend inline bool operator!=(const LocalEntityIterator &lhs,
                                const LocalEntityIterator &rhs) {
    return !(lhs == rhs);
  }

private:
  View _view;
};
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cpppetsc {

template <class MeshType, template <class> class ViewTemplate>
LocalEntityIterator<MeshType, ViewTemplate>::LocalEntityIterator(
    const MeshType *const mesh,
    const typename mesh_type::size_type elementPointIndex)
    : _view(mesh, elementPointIndex) {}

template <class MeshType, template <class> class ViewTemplate>
LocalEntityIterator<MeshType, ViewTemplate> &
LocalEntityIterator<MeshType, ViewTemplate>::operator++() {
  ++_view._entityIndex;
  return *this;
}

template <class MeshType, template <class> class ViewTemplate>
LocalEntityIterator<MeshType, ViewTemplate>
LocalEntityIterator<MeshType, ViewTemplate>::operator++(int) {
  auto copy = *this;
  ++(*this);
  return copy;
}

template <class MeshType, template <class> class ViewTemplate>
LocalEntityIterator<MeshType, ViewTemplate> &
LocalEntityIterator<MeshType, ViewTemplate>::operator--() {
  --_view._entityIndex;
  return *this;
}

template <class MeshType, template <class> class ViewTemplate>
LocalEntityIterator<MeshType, ViewTemplate>
LocalEntityIterator<MeshType, ViewTemplate>::operator--(int) {
  auto copy = *this;
  --(*this);
  return copy;
}

template <class MeshType, template <class> class ViewTemplate>
typename LocalEntityIterator<MeshType, ViewTemplate>::reference
    LocalEntityIterator<MeshType, ViewTemplate>::operator*() const {
  return _view;
}

template <class MeshType, template <class> class ViewTemplate>
typename LocalEntityIterator<MeshType, ViewTemplate>::pointer
    LocalEntityIterator<MeshType, ViewTemplate>::operator->() const {
  return &_view;
}
} // namespace cpppetsc
} // namespace ae108