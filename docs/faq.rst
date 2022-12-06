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

Frequently Asked Questions
==========================

1. How do I limit the number of cores that are used to build the project
   in VS Code?

..

   The `CMake Tools
   extension <https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools>`__
   supports a `number of parameters to configure the
   build <https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/cmake-settings.md>`__.
   Those can be set in the file ``.vscode/settings.json``. To limit the
   number of cores use the parameter ``cmake.parallelJobs``.

2. How do I configure SSH on Windows?

..

   The SSH configuration file on Windows 10 is located at
   ``C:\Users\<username>\.ssh\config``. Note that the contents of this
   file also have an impact on how Docker connects to a remote
   container. For instance, the file could look like this (replace
   ``(...)`` by suitable settings); see also the `Euler
   documentation <https://scicomp.ethz.ch/wiki/Getting_started_with_clusters#SSH>`__:

::

   Host (...)
       User (...)
       HostName (...)
       IdentityFile (...)

..

   In case VS Code asks for the key passphrase too often, consider
   setting up the SSH agent; see the `VS Code
   documentation <https://code.visualstudio.com/docs/remote/troubleshooting#_setting-up-the-ssh-agent>`__
   for more information.

3. I cannot attach to containers via SSH with VS Code.

..

   First check that attaching to the workstation via VS Code and SSH
   works, and that the containers can be started via the
   ``docker-compose.yml`` file in the repository. If this is the case,
   make sure that the rootless Docker context is activated on the
   workstation. For this purpose run ``docker context use rootless`` on
   the workstation. If the problem persists, try uninstalling and
   installing rootless Docker on the workstation (following the
   instructions
   `here <https://docs.docker.com/engine/security/rootless/>`__), and
   then activating the rootless context.

4. VS Code reports a “Bad CMake executable”.

..

   If VS Code shows the error
   ``Bad CMake executable "". Is it installed or settings contain the correct path (cmake.cmakePath)?``,
   then the most likely reason for the error is that the VS Code windows
   is not yet attached to a container. If you are using a custom
   container, then make sure that CMake is installed in the container.

5. VS Code asks for a password even though I have created an SSH key.

..

   Please make sure that you have provided the public key in the
   workstation’s ``authorized_keys`` file; see e.g.
   `scicomp.ethz.ch <https://scicomp.ethz.ch/wiki/Getting_started_with_clusters#SSH_keys>`__.

6. When interacting with ``git`` I am asked for a password even though
   I’ve added a public key.

..

   ``git`` may not be aware of the key to use for the remote repository.
   To solve the problem, add an entry to ``.ssh/config`` (or similar),
   specifying the private key to use (replace ``(...)`` by a suitable
   path):

::

   Host gitlab.ethz.ch
       HostName gitlab.ethz.ch
       IdentityFile (...)

7. Attaching to the container with VS Code fails with the error message
   ``Cannot connect to the Docker daemon at unix:///run/user//docker.sock. Is the docker daemon running?``
   (note the missing user ID between ``/user/`` and ``/docker.sock``).

..

   This error was caused by an incompatibility of the
   `Remote-Containers <https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers>`__
   plugin with the Leibniz configuration. The error should no longer
   occur with the latest version of the plugin.

8. How can I link against a version of PETSc that is compiled with
   enabled debug checks?

..

   Install the Ubuntu package ``libpetsc-real3.12-dbg`` (if you are
   using a real scalar type). Then, before running ``CMake``, set the
   environment variable ``PKG_CONFIG_PATH`` to
   ``/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real-debug/lib/pkgconfig``
   (remember to clear ``CMake's`` cache if it does exist already). If
   you are using VS Code, this can also be accomplished by adding the
   following to ``settings.json``:

::

   "cmake.configureEnvironment": {
     "PKG_CONFIG_PATH": "/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real-debug/lib/pkgconfig"
   }

9. Why do I see the error ``Read -1, expected ..., errno =1`` when
   running an executable in parallel in the Docker container?

..

   This is a known issue, see the `corresponding issue on
   GitHub <https://github.com/open-mpi/ompi/issues/4948>`__. Setting the
   following environment variable works:

::

   export OMPI_MCA_btl_vader_single_copy_mechanism=none

10. How can I add a new developer user to the project?

..

   Please ask the new developer to login to ``gitlab.ethz.ch`` once.
   Then send an email to ``manuel.weberndorfer@id.ethz.ch`` with the
   username of the new developer.

11. The solver does not converge for my scenario. What can I do?

..

   `Jed Brown <https://www.colorado.edu/cs/jed-brown>`__ (one of the
   PETSc developers) offers advice for both `nonlinear
   systems <https://scicomp.stackexchange.com/questions/30/why-is-newtons-method-not-converging>`__
   and `linear
   systems <https://scicomp.stackexchange.com/questions/513/why-is-my-iterative-linear-solver-not-converging>`__.

12. How can I use distributed SuperLU with PETSc?

..

   According to the `PETSc
   documentation <https://petsc.org/release/docs/manualpages/Mat/MATSOLVERSUPERLU_DIST.html#MATSOLVERSUPERLU_DIST>`__,
   this is possible with the flags
   ``-pc_type lu -pc_factor_mat_solver_type superlu_dist``.

13. When I submit my job on Euler via ``bsub``, the error message is
    “Permission denied”. How can I fix this?

..

   Very often, the reason for this error is that the file that ``bsub``
   tries to execute is not executable (e.g. a text file). In this case
   Linux refuses to execute the file. You can confirm that this is the
   issue by attempting to execute the file on the login node
   (i.e. without ``bsub``); this should yield the same error.
   Frequently, the problem is that a different executable was intended
   to be used. Otherwise it is possible to mark a file as executable via
   ``chmod +x`` followed by the path to the file.

14. Docker Desktop fails to start on Windows.

..

   To use Docker Desktop on Windows, make sure to install `Ubuntu 22.04
   from the Microsoft
   store <https://apps.microsoft.com/store/detail/ubuntu-2204-lts/9PN20MSR04DW>`__,
   and activate “Use the WSL 2 based engine” in the Docker settings
   (“General”).
