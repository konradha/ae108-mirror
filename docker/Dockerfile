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

ARG BASE_IMAGE=ubuntu:focal
FROM $BASE_IMAGE

ARG PETSC_SCALAR_TYPE=real

RUN apt-get update && \
        DEBIAN_FRONTEND="noninteractive" \
        TZ="Europe/Zurich" \
        apt-get install -y \
        clang-format \
        cmake \
        g++ \
        git \
        googletest \
        libboost-program-options-dev \
        libeigen3-dev \
        librange-v3-dev \
        ninja-build \
        python3 python3-pip \
        python3-h5py \
        gdb

RUN pip3 install \
        mypy==0.790 \
        black==20.8b1 \
        pylint==2.6.0 \
        sympy==1.7.1 \
        gcovr

RUN apt-get update && \
        DEBIAN_FRONTEND="noninteractive" \
        apt-get install -y \
        libpetsc-$PETSC_SCALAR_TYPE-dev