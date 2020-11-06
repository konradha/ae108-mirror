# © 2020 ETH Zurich, Mechanics and Materials Lab
# © 2020 California Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

ARG BASE_IMAGE=debian:buster
ARG DEBIAN_FRONTEND=noninteractive

FROM $BASE_IMAGE

RUN apt-get update && apt-get install -y clang-format
RUN apt-get update && apt-get install -y cmake
RUN apt-get update && apt-get install -y g++
RUN apt-get update && apt-get install -y googletest
RUN apt-get update && apt-get install -y libboost-program-options-dev
RUN apt-get update && apt-get install -y libeigen3-dev
RUN apt-get update && apt-get install -y ninja-build
RUN apt-get update && apt-get install -y petsc-dev
RUN apt-get update && apt-get install -y python3 python3-pip

RUN pip3 install gcovr