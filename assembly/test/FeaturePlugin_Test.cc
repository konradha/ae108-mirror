// © 2022 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/assembly/FeaturePlugin.h"

namespace {

DEFINE_ASSEMBLER_PLUGIN(TestPlugin, testMethod, (double)) {}
DEFINE_CONST_ASSEMBLER_PLUGIN(TestPluginConst, testMethodConst, (double)) {}

} // namespace

namespace ae108 {
namespace assembly {
namespace {

DEFINE_ASSEMBLER_PLUGIN(TestPlugin, testMethod, (double)) {}
DEFINE_CONST_ASSEMBLER_PLUGIN(TestPluginConst, testMethodConst, (double)) {}

} // namespace
} // namespace assembly
} // namespace ae108
