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

ARG BASE_IMAGE=ubuntu:jammy
FROM $BASE_IMAGE

ARG PETSC_SCALAR_TYPE=real

RUN apt-get update && \
  DEBIAN_FRONTEND="noninteractive" \
  TZ="Europe/Zurich" \
  apt-get install -y \
  clang-format \
  cmake \
  g++ \
  gdb \
  git \
  libboost-program-options-dev \
  libeigen3-dev \
  libgmock-dev \
  libgtest-dev \
  librange-v3-dev \
  libxrender1 \
  ninja-build \
  python3 python3-pip \
  python3-scipy

RUN pip3 install \
  mypy==0.931 \
  black==22.3.0 \
  pylint==2.12.2 \
  sympy==1.7.1 \
  gcovr==5.0 \
  h5py==3.6.0 \
  jsonschema==4.6.0

RUN pip3 install vtk --extra-index-url https://gitlab.kitware.com/api/v4/projects/13/packages/pypi/simple

RUN apt-get update && \
  DEBIAN_FRONTEND="noninteractive" \
  apt-get install -y \
  libpetsc-${PETSC_SCALAR_TYPE}-dev \
  libslepc-${PETSC_SCALAR_TYPE}3.15-dev
