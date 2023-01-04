# © 2020 ETH Zurich, Mechanics and Materials Lab
# © 2020 California Institute of Technology
#
# This file is part of ae108.
#
# ae108 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any
# later version.
#
# ae108 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ae108. If not, see <https://www.gnu.org/licenses/>.

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
  jsonschema==4.6.0 \
  vtk==9.2.2

RUN apt-get update && \
  DEBIAN_FRONTEND="noninteractive" \
  apt-get install -y \
  libpetsc-${PETSC_SCALAR_TYPE}-dev \
  libslepc-${PETSC_SCALAR_TYPE}3.15-dev
