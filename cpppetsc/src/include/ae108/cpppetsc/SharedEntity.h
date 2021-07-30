// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include <memory>
#include <petscsys.h>
#include <type_traits>

namespace ae108 {
namespace cpppetsc {

template <typename T, typename Policy> class SharedEntity {
public:
  /**
   * @brief Takes ownership of `t` without increasing the reference count.
   */
  explicit SharedEntity(T t);

  /**
   * @brief Creates a shallow copy of the entity by increasing the reference
   * count.
   */
  SharedEntity(const SharedEntity &in);

  /**
   * @brief Creates a shallow copy of the entity by increasing the reference
   * count.
   */
  SharedEntity &operator=(const SharedEntity &in);

  /**
   * @brief Takes ownership without increasing the reference count.
   */
  SharedEntity(SharedEntity &&);

  /**
   * @brief Takes ownership without increasing the reference count.
   */
  SharedEntity &operator=(SharedEntity &&);

  /**
   * @brief Returns the entity.
   */
  T get() const noexcept;

  /**
   * @brief Returns true if and only if the contained entity is not the nullptr.
   */
  explicit operator bool() const noexcept;

  using size_type = PetscInt;

  /**
   * @brief Returns the reference count of the entity.
   */
  size_type count() const;

private:
  /**
   * @brief Decreases the reference count.
   */
  static void destroy(T t) noexcept;

  static_assert(std::is_pointer<T>::value, "Only pointers T are supported.");
  std::unique_ptr<typename std::remove_pointer<T>::type, decltype(&destroy)> t_;
};

enum class OwnershipType { Take, Share };

/**
 * @brief Create a SharedEntity.
 *
 * @param own Take: do not increase reference count;
 * Share: increase reference count
 */
template <typename Policy, typename T>
SharedEntity<T, Policy>
makeSharedEntity(T t, const OwnershipType own = OwnershipType::Take);

} // namespace cpppetsc
} // namespace ae108

namespace ae108 {
namespace cpppetsc {

template <class T, class Policy>
SharedEntity<T, Policy>::SharedEntity(T t) : t_(t, &destroy) {}

template <class T, class Policy>
SharedEntity<T, Policy>::SharedEntity(const SharedEntity &in)
    : t_((Policy::handleError(PetscObjectReference((PetscObject)(in.get()))),
          in.get()),
         &destroy) {}

template <class T, class Policy>
SharedEntity<T, Policy> &
SharedEntity<T, Policy>::operator=(const SharedEntity &in) {
  auto copy = in;
  using std::swap;
  swap(*this, copy);
  return *this;
}

template <class T, class Policy>
SharedEntity<T, Policy>::SharedEntity(SharedEntity &&) = default;

template <class T, class Policy>
SharedEntity<T, Policy> &
SharedEntity<T, Policy>::operator=(SharedEntity &&) = default;

template <class T, class Policy>
SharedEntity<T, Policy>::operator bool() const noexcept {
  return t_.get();
}

template <class T, class Policy>
T SharedEntity<T, Policy>::get() const noexcept {
  return t_.get();
}

template <class T, class Policy>
typename SharedEntity<T, Policy>::size_type
SharedEntity<T, Policy>::count() const {
  auto result = size_type{0};
  Policy::handleError(PetscObjectGetReference((PetscObject)t_.get(), &result));
  return result;
}

template <class T, class Policy>
void SharedEntity<T, Policy>::destroy(T t) noexcept {
  Policy::handleError(PetscObjectDereference((PetscObject)t));
}

template <typename Policy, typename T>
SharedEntity<T, Policy> makeSharedEntity(T t, const OwnershipType own) {
  return SharedEntity<T, Policy>{
      own == OwnershipType::Take
          ? t
          : (Policy::handleError(PetscObjectReference((PetscObject)t)), t)};
}

} // namespace cpppetsc
} // namespace ae108