.. raw:: html

   <!---
    © 2022 ETH Zurich, Mechanics and Materials Lab

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
   -->

Quick start with Docker
=======================

There are different options to deploy this framework. For a quick start
and versatile development environment that works on many platforms –
from your local laptop to large scale high performance computers – we
rely on Docker.

Why Docker?
-----------

Configuring a workstation directly and working on this machine has
disadvantages that become problematic over time.

-  **Lock-in**: The exact configuration cannot be reproduced (easily),
   so there’s no other option than to keep developing on that single
   workstation, even if there are better options (e.g. a second faster
   or less occupied workstation).
-  **Reproducibility**: Since the configuration is not documented,
   reproducing a simulation is hard if the configuration has been
   changed.
-  **Updates**: It is difficult to update the environment
   (e.g. libraries), since older versions of the code may depend on
   older versions of the library. This leads to broken builds without
   changed code for older versions of the code.
-  **Permissions**: Installing libraries, even if they are available via
   ``apt-get``, requires support by an administrator. It is difficult to
   quickly test a library.
-  **Testability**: Automatic tools that check the code also need to run
   on this workstation to verify that it works.

This leads to a situation where it’s hard to improve both the
environment and the code.

The benefit of using Docker is to have a development environment that
can be shared between developers, can be created on all platforms, and
can be modified easily. The commands that define the environment are in
a single short file, the ``Dockerfile``. This Dockerfile is `versioned
together with the
code <https://gitlab.ethz.ch/mechanics-and-materials/ae108/-/blob/master/docker/Dockerfile>`__.

Requirements
------------

Docker supports Linux, MacOS, and Windows hosts. In the following, one
of the following operating systems is required:

-  Ubuntu 20.04 Focal
-  Windows 10/11
-  MacOS 12 Monterey

Tutorial
--------

The following tutorial guides you through a quick way to set up a
development environment on your machine using
`Docker <https://www.docker.com/>`__ and `Visual Studio
Code <https://code.visualstudio.com/>`__.

1. Install Git by following the instructions available at
   `git-scm.com <https://git-scm.com/downloads>`__ (use the default
   settings on Windows).
2. Install Docker by following the instructions available at
   `docker.com <https://docs.docker.com/get-docker/>`__.

   If you are using Linux then install rootless Docker via the
   `instructions <https://docs.docker.com/engine/security/rootless/>`__
   and activate it using:

   .. code:: bash

      docker context use rootless

3. If you are on Linux, install Docker Compose by following the
   instructions available at
   `docker.com <https://docs.docker.com/compose/install/>`__ (Docker
   Desktop for Windows/Mac already comes with Docker Compose). In
   particular, the package ``docker-compose`` can also be installed via
   ``pip``:

   .. code:: bash

      pip3 install docker-compose

4. Install Visual Studio Code by following the instructions available at
   `visualstudio.com <https://code.visualstudio.com/>`__. In addition,
   install the following Visual Studio code extensions (by Microsoft):

   -  `C/C++ <https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools>`__
   -  `Docker <https://marketplace.visualstudio.com/items?itemName=ms-azuretools.vscode-docker>`__
   -  `Remote -
      Containers <https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers>`__
   -  `CMake
      Tools <https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools>`__

5. Click the “Clone” button on the `GitLab
   website <https://gitlab.ethz.ch/mechanics-and-materials/ae108>`__ of
   the repository and select “Visual Studio Code” in the menu. After
   you’ve chosen a target directory a new Visual Studio Code window
   should open automatically. Confirm opening the new folder when VS
   Code asks you to.
6. VS Code shows a dialog box informing you that the repository provides
   a “Dev Container” configuration. Start this container via the button
   “Reopen in Container”. This can take a few minutes to complete.
7. Click the CMake icon in the menu on the left, and click the
   “Configure all Projects” button on the top of the CMake pane.
8. Build the project by clicking “Build all Projects” on the top of the
   CMake pane.

Execute an interactive bash shell in the container
--------------------------------------------------

After starting the container as described above, you can also jump right
into the container from your terminal (outside of VS Code) by executing
an interactive bash shell. The command is:

.. code:: bash

   docker exec -it <container name> bash

If you do not know your ``<container name>``, you may always look it up
with the following command:

.. code:: bash

   docker ps
