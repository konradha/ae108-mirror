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
