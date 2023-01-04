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

#include "ae108/assembly/AssemblerTypeTraits.h"
#include <functional>
#include <tuple>
#include <utility>

#define DEFINE_ASSEMBLER_METHOD_BASE(PluginName, methodName, cvQualifiers)     \
  template <class... Args> void methodName(Args &&...args) cvQualifiers {      \
    if constexpr (!::ae108::assembly::IsGroupTypeTrait<                        \
                      typename PluginName::assembler_type>::value) {           \
      execute(std::forward<Args>(args)...);                                    \
    } else {                                                                   \
      const auto conditionalCall = [&](auto &a) {                              \
        using A = std::decay_t<decltype(a)>;                                   \
        if constexpr (std::is_base_of_v<PluginName<A>, A>) {                   \
          static_cast<cvQualifiers PluginName<A> &>(a).execute(args...);       \
        }                                                                      \
      };                                                                       \
      std::apply(                                                              \
          [&](auto &...as) { (..., std::invoke(conditionalCall, as)); },       \
          this->assembler().assemblers());                                     \
    }                                                                          \
  }

#define DEFINE_ASSEMBLER_PLUGIN_BASE(PluginName, methodName, methodArguments,  \
                                     cvQualifiers)                             \
  template <class Assembler>                                                   \
  class PluginName                                                             \
      : public ::ae108::assembly::FeaturePlugin<Assembler, PluginName> {       \
    using assembler_type = Assembler;                                          \
    using element_type =                                                       \
        typename ::ae108::assembly::ElementTypeTrait<Assembler>::type;         \
    using matrix_type =                                                        \
        typename ::ae108::assembly::MatrixTypeTrait<Assembler>::type;          \
    using mesh_type =                                                          \
        typename ::ae108::assembly::MeshTypeTrait<Assembler>::type;            \
    using size_type =                                                          \
        typename ::ae108::assembly::SizeTypeTrait<Assembler>::type;            \
    using value_type =                                                         \
        typename ::ae108::assembly::ValueTypeTrait<Assembler>::type;           \
    using real_type =                                                          \
        typename ::ae108::assembly::RealTypeTrait<Assembler>::type;            \
    using vector_type =                                                        \
        typename ::ae108::assembly::VectorTypeTrait<Assembler>::type;          \
                                                                               \
  public:                                                                      \
    void execute methodArguments cvQualifiers;                                 \
    DEFINE_ASSEMBLER_METHOD_BASE(PluginName, methodName, cvQualifiers)         \
                                                                               \
  protected:                                                                   \
    ~PluginName() = default;                                                   \
  };                                                                           \
                                                                               \
  template <class Assembler>                                                   \
  void PluginName<Assembler>::execute methodArguments cvQualifiers

#define DEFINE_CONST_ASSEMBLER_PLUGIN(PluginName, methodName, methodArguments) \
  DEFINE_ASSEMBLER_PLUGIN_BASE(PluginName, methodName, methodArguments, const)

#define DEFINE_ASSEMBLER_PLUGIN(PluginName, methodName, methodArguments)       \
  DEFINE_ASSEMBLER_PLUGIN_BASE(PluginName, methodName, methodArguments, )

namespace ae108 {
namespace assembly {

template <class Assembler_, template <class> class Plugin_>
class FeaturePlugin {
  using assembler_type = Assembler_;
  template <class T> using plugin_type = Plugin_<T>;

protected:
  /**
   * @brief Returns a reference to the Assembler base class.
   */
  assembler_type &assembler();

  /**
   * @brief Returns a reference to the Assembler base class.
   */
  const assembler_type &assembler() const;

  ~FeaturePlugin() = default;
};
} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace assembly {

template <class Assembler, template <class> class Plugin>
Assembler &FeaturePlugin<Assembler, Plugin>::assembler() {
  return static_cast<Assembler &>(static_cast<Plugin<Assembler> &>(*this));
}

template <class Assembler, template <class> class Plugin>
const Assembler &FeaturePlugin<Assembler, Plugin>::assembler() const {
  return static_cast<const Assembler &>(
      static_cast<const Plugin<Assembler> &>(*this));
}
} // namespace assembly
} // namespace ae108
