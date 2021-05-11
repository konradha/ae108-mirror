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

#include <cassert>
#include <cstddef>
#include <iterator>
#include <type_traits>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief A bidirectional iterator to index-access a range.
 */
template <class RangeType, class IndexType = typename RangeType::size_type,
          class ValueType = typename RangeType::value_type>
class IndexBasedAccessIterator {
  static_assert(std::is_integral<IndexType>::value,
                "IndexType is an integral type.");

public:
  using value_type = ValueType;
  using iterator_category = std::bidirectional_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using reference = value_type;
  struct PointerProxy;
  using pointer = PointerProxy;

  IndexBasedAccessIterator() = default;
  explicit IndexBasedAccessIterator(const RangeType *const range,
                                    const IndexType index);

  IndexBasedAccessIterator &operator++();

  IndexBasedAccessIterator operator++(int);

  IndexBasedAccessIterator &operator--();

  IndexBasedAccessIterator operator--(int);

  reference operator*() const;
  pointer operator->() const;

  friend inline bool operator==(const IndexBasedAccessIterator &lhs,
                                const IndexBasedAccessIterator &rhs) {
    return lhs._range == rhs._range && lhs._index == rhs._index;
  }

  friend inline bool operator!=(const IndexBasedAccessIterator &lhs,
                                const IndexBasedAccessIterator &rhs) {
    return !(lhs == rhs);
  }

private:
  const RangeType *_range = nullptr;
  IndexType _index = 0;
};
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cpppetsc {

template <class RangeType, class IndexType, class ValueType>
IndexBasedAccessIterator<RangeType, IndexType, ValueType>::
    IndexBasedAccessIterator(const RangeType *const range,
                             const IndexType index)
    : _range(range), _index(index) {}

template <class RangeType, class IndexType, class ValueType>
IndexBasedAccessIterator<RangeType, IndexType, ValueType> &
IndexBasedAccessIterator<RangeType, IndexType, ValueType>::operator++() {
  ++_index;
  return *this;
}

template <class RangeType, class IndexType, class ValueType>
IndexBasedAccessIterator<RangeType, IndexType, ValueType>
IndexBasedAccessIterator<RangeType, IndexType, ValueType>::operator++(int) {
  auto copy = *this;
  ++(*this);
  return copy;
}

template <class RangeType, class IndexType, class ValueType>
IndexBasedAccessIterator<RangeType, IndexType, ValueType> &
IndexBasedAccessIterator<RangeType, IndexType, ValueType>::operator--() {
  --_index;
  return *this;
}

template <class RangeType, class IndexType, class ValueType>
IndexBasedAccessIterator<RangeType, IndexType, ValueType>
IndexBasedAccessIterator<RangeType, IndexType, ValueType>::operator--(int) {
  auto copy = *this;
  --(*this);
  return copy;
}

template <class RangeType, class IndexType, class ValueType>
typename IndexBasedAccessIterator<RangeType, IndexType, ValueType>::reference
IndexBasedAccessIterator<RangeType, IndexType, ValueType>::operator*() const {
  assert(_range);
  return (*_range)[_index];
}

template <class RangeType, class IndexType, class ValueType>
struct IndexBasedAccessIterator<RangeType, IndexType, ValueType>::PointerProxy {
  ValueType value;
  const ValueType *operator->() { return &value; }
};

template <class RangeType, class IndexType, class ValueType>
typename IndexBasedAccessIterator<RangeType, IndexType, ValueType>::pointer
IndexBasedAccessIterator<RangeType, IndexType, ValueType>::operator->() const {
  return PointerProxy{*(*this)};
}
} // namespace cpppetsc
} // namespace ae108