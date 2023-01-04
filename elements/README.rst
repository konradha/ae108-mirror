..
    © 2020 ETH Zurich, Mechanics and Materials Lab

    This file is part of ae108.

    ae108 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any
    later version.

    ae108 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ae108. If not, see <https://www.gnu.org/licenses/>.

The Elements Library
====================

Introduction
------------

This library provides a framework to write elements and material models together with examples.
For instance, the library contains an implementation of the ``Hexa8`` element with a Hookean material model.

The main idea behind this library is that an element is just a material model (a way to compute an entity like energy) combined with an integrator (a way to integrate this entity over the element).
These two concepts are (almost) orthogonal and the two components can therefore be chosen (almost) independently.

Components
----------

Currently this library has eight parts.

Elements
~~~~~~~~

The files for this part are in the base ``src`` folder.
It contains traits to compute energy, forces, and stiffness matrix of an element.
There are also traits to compute the forces and the stiffness matrix numerically.

These traits are implemented for the ``CoreElement`` element (which should be applicable to many use cases) and the ``Minimal`` element (which shows how to implement a valid element with minimal effort).

Material Models
~~~~~~~~~~~~~~~

This subfolder contains the traits to compute energy, stress, strain, and tangent matrix of a material model.
There are also traits to compute the stress and the tangent matrix numerically.

These traits are implemented for the ``Hookean`` material model and the ``Minimal`` material model (which shows how to implement a valid material model with minimal effort).

Integrator
~~~~~~~~~~

This subfolder contains a trait to integrate an entity over an element.

There’s also an implementation of an ``IsoparametricIntegrator``, which uses an isoparametric embedding, a set of shape functions, and a quadrature rule in reference space to integrate entities.

Embedding
~~~~~~~~~

This subfolder contains a trait to embed a point in reference space into physical space and a trait to compute the Jacobian of this transformation.
There is also a trait to compute the Jacobian numerically.

These traits are implemented for the ``IsoparametricEmbedding`` which uses a set of shape functions to embed an element into physical space.

Shape
~~~~~

This subfolder contains traits to compute the values, gradients, and support points of a set of shape functions.
Here, “support point” refers to a point where a shape function is equal to one whereas all other shape functions are equal to zero.
There is also a trait to compute the gradients numerically.

The part includes three sets of shape functions for which these traits are implemented: ``Hexa8``, ``Quad4``, and ``Seg2``.

Quadrature
~~~~~~~~~~

This subfolder contains a trait to integrate a function over a reference element.

This trait is implemented for a set of quadrature rules and a script to generate them automatically.

Tensor
~~~~~~

This subfolder contains utilities for working with C++ arrays. For instance, it contains functions that interpret such arrays as Eigen matrices.

Mesh
~~~~

This subfolder contains functions that generate meshes.

FAQ
---

What do I need to do to write a material model?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a header file in ``src/include/ae108/elements/materialmodels/`` and an implementation file (that includes the header file) in ``src/materialmodels/`` and add the latter to the ``CMakeLists.txt`` in the ``src`` directory.
If you don’t want your material model to be part of this library you can also put those files at a different location.

Now write a class that contains all the data and implement the traits that you want to support.
Note that the new material model should derive from ``MaterialModelBase``.
Feel free to look at the ``Hookean`` material model for inspiration.

To test your material model add an implementation file to the ``test/materialmodels`` directory and add it to the ``CMakeLists.txt`` int the ``test`` directory.
Make sure you instantiate the ``MaterialModel_Test`` set of tests for your class to run many checks that apply to all material models (e.g. consistent derivatives) automatically.
Again, have a look at ``Hookean_Test`` to see an example of how this can be done.

Why are there traits (like ``ComputeEnergyTrait``) and free functions (like ``compute_energy``)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Traits make it very easy to choose an implementation of “computing energy” without copying and pasting code.
This provides a lot of flexibility when writing elements.
However, this flexibility may make it harder to use an element because every specialization of an ``ComputeEnergyTrait`` could e.g. return a different type.
So, the goal is to only depend on properties of ``ElementBase`` so that calling code can adjust easily.

To solve this problem the free function ``compute_energy`` promises an interface to users: it always returns ``ElementBase::Energy`` which is a name for the element’s ``value_type``.
It also simplifies using the trait.

To sum up:

-  Authors of elements should implement traits that are compatible with the free functions.
-  Users of elements should use the free functions to interact with elements.

The same reasoning applies to other components like material models, etc.

Why are the ``operator()`` methods of many traits template functions?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The idea behind using traits is to make it possible to reuse implemented traits for different classes.
If an implementer of e.g. a material model ``B`` decides to reuse the implementation of the ``ComputeEnergyTrait`` of class ``A``, the following code should be sufficient:

.. code-block:: cpp

   template<>
   struct ComputeEnergyTrait<B> : ComputeEnergyTrait<A> {};

Now assume that you have an instance ``b`` of a ``B`` and that ``operator()`` is not a template.
When calling the free function ``compute_energy(b, ...)`` compilation would fail because ``operator()`` of ``B``'s ``ComputeEnergyTrait`` needs to be called with an instance of ``A`` (since we are inheriting this operator from ``A``'s trait).

By declaring ``operator()`` a template function this use case works.
