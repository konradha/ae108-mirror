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

#include "ae108/assembly/AssemblerTypeTraits.h"
#include <utility>

#define DEFINE_ASSEMBLER_METHOD_BASE(methodName, cvQualifiers)                 \
  template <class... Args> void methodName(Args &&... args) cvQualifiers {     \
    this->dispatch(std::forward<Args>(args)...);                               \
  }

#define DEFINE_ASSEMBLER_PLUGIN_BASE(PluginName, methodName, methodArguments,  \
                                     cvQualifiers)                             \
  template <class Assembler>                                                   \
  class PluginName                                                             \
      : public ::ae108::assembly::FeaturePlugin<Assembler, PluginName> {       \
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
    DEFINE_ASSEMBLER_METHOD_BASE(methodName, cvQualifiers)                     \
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

template <class Assembler, template <class> class Plugin> class FeaturePlugin {
protected:
  /**
   * @brief Returns a reference to the Assembler base class.
   */
  Assembler &assembler();

  /**
   * @brief Returns a reference to the Assembler base class.
   */
  const Assembler &assembler() const;

  /**
   * @brief Calls execute with the given parameters for all (member) assemblers
   * that have the plugin (const version).
   */
  template <class... Args> void dispatch(Args &&... args) const;

  /**
   * @brief Calls execute with the given parameters for all (member) assemblers
   * that have the plugin (nonconst version).
   */
  template <class... Args> void dispatch(Args &&... args);

  ~FeaturePlugin() = default;
};
} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include "ae108/assembly/utilities/StaticLooper.h"
#include <type_traits>

namespace ae108 {
namespace assembly {

namespace detail {
template <class Assembler, template <class> class Plugin>
using HasPlugin = std::is_base_of<Plugin<Assembler>, Assembler>;

template <class Assembler, template <class> class Plugin>
using PluginType = Plugin<typename std::remove_const<Assembler>::type>;

template <class Assembler, template <class> class Plugin>
using CastType =
    typename std::conditional<std::is_const<Assembler>::value,
                              const PluginType<Assembler, Plugin> &,
                              PluginType<Assembler, Plugin> &>::type;

template <template <class> class Plugin> struct Executor {
  template <class Assembler, class... Args>
  void operator()(Assembler &assembler, Args &&... args) const {
    operatorImpl(
        HasPlugin<typename std::remove_const<Assembler>::type, Plugin>{},
        assembler, std::forward<Args>(args)...);
  }

private:
  template <class Assembler, class... Args>
  void operatorImpl(std::true_type, Assembler &assembler,
                    Args &&... args) const {
    static_cast<CastType<Assembler, Plugin>>(assembler).execute(
        std::forward<Args>(args)...);
  }

  template <class Assembler, class... Args>
  void operatorImpl(std::false_type, Assembler &, Args &&...) const {}
};

/**
 * @brief Implementation of dispatch for the assembler-group case.
 */
template <template <class> class Plugin, class Assembler, class... Args>
void dispatchImpl(std::true_type, Assembler &assembler, Args &&... args) {
  utilities::StaticLooper<NumberOfMembersTypeTrait<typename std::remove_const<
      Assembler>::type>::value>::run(detail::Executor<Plugin>{}, assembler,
                                     std::forward<Args>(args)...);
}

/**
 * @brief Implementation of dispatch for the non-assembler-group case.
 */
template <template <class> class Plugin, class Assembler, class... Args>
void dispatchImpl(std::false_type, Assembler &assembler, Args &&... args) {
  static_cast<CastType<Assembler, Plugin>>(assembler).execute(
      std::forward<Args>(args)...);
}
} // namespace detail

template <class Assembler, template <class> class Plugin>
Assembler &FeaturePlugin<Assembler, Plugin>::assembler() {
  return static_cast<Assembler &>(static_cast<Plugin<Assembler> &>(*this));
}

template <class Assembler, template <class> class Plugin>
const Assembler &FeaturePlugin<Assembler, Plugin>::assembler() const {
  return static_cast<const Assembler &>(
      static_cast<const Plugin<Assembler> &>(*this));
}

template <class Assembler, template <class> class Plugin>
template <class... Args>
void FeaturePlugin<Assembler, Plugin>::dispatch(Args &&... args) const {
  detail::dispatchImpl<Plugin>(IsGroupTypeTrait<Assembler>{}, assembler(),
                               std::forward<Args>(args)...);
}

template <class Assembler, template <class> class Plugin>
template <class... Args>
void FeaturePlugin<Assembler, Plugin>::dispatch(Args &&... args) {
  detail::dispatchImpl<Plugin>(IsGroupTypeTrait<Assembler>{}, assembler(),
                               std::forward<Args>(args)...);
}
} // namespace assembly
} // namespace ae108
