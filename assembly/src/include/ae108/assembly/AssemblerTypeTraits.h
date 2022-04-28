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

#include <cstddef>
#include <type_traits>

namespace ae108 {
namespace assembly {

/**
 * @brief Deduces the policy type corresponding to an Assembler class.
 */
template <class Assembler> struct PolicyTypeTrait {};

/**
 * @brief Deduces the mesh type corresponding to an Assembler class.
 */
template <class Assembler> struct MeshTypeTrait {};

/**
 * @brief Deduces the vector type corresponding to an Assembler class.
 */
template <class Assembler> struct VectorTypeTrait {};

/**
 * @brief Deduces the matrix type corresponding to an Assembler class.
 */
template <class Assembler> struct MatrixTypeTrait {};

/**
 * @brief Deduces the size type corresponding to an Assembler class.
 */
template <class Assembler> struct SizeTypeTrait {};

/**
 * @brief Deduces the value type corresponding to an Assembler class.
 */
template <class Assembler> struct ValueTypeTrait {};

/**
 * @brief Deduces the real type corresponding to an Assembler class.
 */
template <class Assembler> struct RealTypeTrait {};

/**
 * @brief Deduces the element type corresponding to an Assembler class.
 */
template <class Assembler> struct ElementTypeTrait {};

/**
 * @brief Deduces the plugin type corresponding to an Assembler class.
 */
template <class Assembler> struct PluginTypeTrait {};

/**
 * @brief Deduces whether the assembler is an assembler group. By default
 * (designed for the nongroup case) the value is false.
 */
template <class Assembler> struct IsGroupTypeTrait : std::false_type {};
} // namespace assembly
} // namespace ae108
