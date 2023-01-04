// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#pragma once

#include <type_traits>
#include <utility>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief A wrapper that adds a tag to a given class.
 */
template <class Entity, class Tag> class TaggedEntity {
public:
  using value_type = typename std::remove_cv<
      typename std::remove_reference<Entity>::type>::type;
  using tag_type = Tag;

  /**
   * @brief Constructs the contained value with the given arguments.
   */
  template <class... Args> explicit constexpr TaggedEntity(Args &&...args);

  /**
   * @brief Calls operator() of the contained value.
   */
  template <class... Args>
  auto operator()(Args &&...args)
      -> decltype(std::declval<value_type>()(std::forward<Args>(args)...));

  /**
   * @brief Calls operator() of the contained value.
   */
  template <class... Args>
  constexpr auto operator()(Args &&...args) const
      -> decltype(std::declval<const value_type>()(
          std::forward<Args>(args)...));

  /**
   * @brief Returns a const reference to the contained value.
   */
  constexpr const value_type &unwrap() const &noexcept;

  /**
   * @brief Returns a reference to the contained value.
   */
  value_type &unwrap() &noexcept;

  /**
   * @brief Returns an r-value reference to the contained value.
   */
  value_type &&unwrap() &&noexcept;

  /**
   * @brief Automatically unwraps the contained value.
   */
  constexpr operator const value_type &() const &noexcept;

  /**
   * @brief Automatically unwraps the contained value.
   */
  operator value_type &() &noexcept;

  /**
   * @brief Automatically unwraps the contained value.
   */
  operator value_type &&() &&noexcept;

private:
  value_type _value;
};

/**
 * @brief Wraps the given entity in a TaggedEntity with the provided tag.
 */
template <class Tag, class Entity>
constexpr TaggedEntity<Entity, Tag> tag(Entity &&entity);
} // namespace cpppetsc
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cpppetsc {

template <class Entity, class Tag>
template <class... Args>
constexpr TaggedEntity<Entity, Tag>::TaggedEntity(Args &&...args)
    : _value(std::forward<Args>(args)...) {}

template <class Entity, class Tag>
template <class... Args>
auto TaggedEntity<Entity, Tag>::operator()(Args &&...args)
    -> decltype(std::declval<typename TaggedEntity<Entity, Tag>::value_type>()(
        std::forward<Args>(args)...)) {
  return _value(std::forward<Args>(args)...);
}

template <class Entity, class Tag>
template <class... Args>
constexpr auto TaggedEntity<Entity, Tag>::operator()(Args &&...args) const
    -> decltype(std::declval<
                const typename TaggedEntity<Entity, Tag>::value_type>()(
        std::forward<Args>(args)...)) {
  return _value(std::forward<Args>(args)...);
}

template <class Entity, class Tag>
constexpr const typename TaggedEntity<Entity, Tag>::value_type &
TaggedEntity<Entity, Tag>::unwrap() const &noexcept {
  return _value;
}

template <class Entity, class Tag>
typename TaggedEntity<Entity, Tag>::value_type &
TaggedEntity<Entity, Tag>::unwrap() &noexcept {
  return _value;
}

template <class Entity, class Tag>
typename TaggedEntity<Entity, Tag>::value_type &&
TaggedEntity<Entity, Tag>::unwrap() &&noexcept {
  return std::move(_value);
}

template <class Entity, class Tag>
constexpr TaggedEntity<Entity, Tag>::operator const typename TaggedEntity<
    Entity, Tag>::value_type &() const &noexcept {
  return _value;
}

template <class Entity, class Tag>
TaggedEntity<Entity, Tag>::operator typename TaggedEntity<
    Entity, Tag>::value_type &() &noexcept {
  return _value;
}

template <class Entity, class Tag>
TaggedEntity<Entity, Tag>::operator typename TaggedEntity<
    Entity, Tag>::value_type &&() &&noexcept {
  return std::move(_value);
}

template <class Tag, class Entity>
constexpr TaggedEntity<Entity, Tag> tag(Entity &&entity) {
  return TaggedEntity<Entity, Tag>(std::forward<Entity>(entity));
}
} // namespace cpppetsc
} // namespace ae108