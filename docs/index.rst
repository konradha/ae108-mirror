..
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

:html_theme.sidebar_secondary.remove:

===============================================
ae108: a scalable variational framework for FEA
===============================================

.. only:: html
  
  .. grid:: 5
    :gutter: 0
    :margin: 0
    :padding: 0

    .. grid-item::

      .. image:: _static/trussIndentation.gif
        :target: examples/README.html
        :class: dark-light p-2

    .. grid-item::

      .. image:: _static/crystalPlasticity.png 
        :target: examples/README.html
        :class: dark-light p-2

    .. grid-item::

      .. image:: _static/wavePropagation.gif 
        :target: examples/README.html
        :class: dark-light p-2

    .. grid-item::

      .. image:: _static/homogenization.png 
        :target: examples/README.html
        :class: dark-light p-2

    .. grid-item::

      .. image:: _static/neoHookeanBall.gif 
        :target: examples/README.html
        :class: dark-light p-2

*ae108*  is a scalable (parallel) C++ framework for computational solid mechanics simulations using a variational approach, primarily focusing on the Finite Element Analysis (FEA).
See the `Computational Solid Mechanics lecture notes <https://ethz.ch/content/dam/ethz/special-interest/mavt/mechanical-systems/mm-dam/documents/Notes/CompMech_Notes.pdf>`_  for more information about the approach that is used.
The code is developed by the `Mechanics & Materials Lab at ETH Zurich <https://mm.ethz.ch/>`_.

.. note::
  *ae108* is a research library. The libary, as well as this documentation, is under active development. The contents are constantly expanding and may be subject to frequent changes. For the latest information, please refer to our `GitLab <https://gitlab.ethz.ch/mechanics-and-materials/ae108>`_.

For general citations on *ae108* please use the following:

.. literalinclude:: _static/ae108.bib
   :language: none
   :start-at: @misc{ae108
   :end-at: year
   :append: }

.. toctree::
  :hidden:

  Quick Start <quickstart.rst>
  Examples <examples/README.rst>
  Libraries <libraries.rst>
  FAQ <faq.rst>

  