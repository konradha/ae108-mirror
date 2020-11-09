..
    © 2020 ETH Zurich, Mechanics and Materials Lab

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

The CmdLine Library
===================

Introduction
------------

This library assists with reading parameters from the command line.
It supports different types of arguments, including floating point numbers, integers, and booleans.
Moreover, it automatically generates a ``--help`` flag that provides custom information about the available flags.

Tutorial
--------

To use this library, you only need to familiarize yourself with a single class.
This class is called ``CommandLineOptionParser``.
It's in the namespace ``ae108::cmdline``.

Example Application
^^^^^^^^^^^^^^^^^^^

Let's try using this class in a simple example.
In this example, we are going greet the world if greetings are enabled via a command line flag.

First, we include the necessary headers and start with defining a ``main`` function.

.. literalinclude:: cmdline-example.cc
    :language: cpp
    :lines: 15-19

To parse the command line parameters, we'll construct a ``CommandLineOptionParser`` by providing a stream.
Let's say we want to print error messages and warnings to ``stderr``.
In addition, we define a boolean variable ``enable_greeting`` that will store whether the greeting was enabled.

.. literalinclude:: cmdline-example.cc
    :language: cpp
    :lines: 21-22

We'll add a flag ``--enable_greeting`` (with a short form ``-g``) together with a help text describing the flag.
The ``CommandLineOptionParser`` class provides a ``withOption`` method that we are going to use to achieve that.

.. literalinclude:: cmdline-example.cc
    :language: cpp
    :lines: 23

Now that we've configured the ``CommandLineOptionParser``, the only thing that is missing is using it to parse the command line options in ``argc``/``argv``.

.. literalinclude:: cmdline-example.cc
    :language: cpp
    :lines: 24-25

Finally we print a message if greetings have been enabled.

.. literalinclude:: cmdline-example.cc
    :language: cpp
    :lines: 26-

Let's see the application in action.

.. code-block:: shell-session

    $ cmdline/examples/ae108-CmdLineExample
    $ cmdline/examples/ae108-CmdLineExample --enable_greeting=yes
    Hello world!
    $ cmdline/examples/ae108-CmdLineExample --enable_greeting=no
    $ cmdline/examples/ae108-CmdLineExample -g yes
    Hello world!

As promised, there's also a ``--help`` flag:

.. code-block:: shell-session

    $ cmdline/examples/ae108-CmdLineExample --help
    Command line options:
    -h [ --help ]                Show this help.
    -g [ --enable_greeting ] arg Print a greeting.

Moreover, a warning message is printed out to ``stderr`` if an unknown flag is used.

.. code-block:: shell-session

    $ cmdline/examples/ae108-CmdLineExample --unknown_flag=123
    Warning: The following options were not recognized by cmdline: '--unknown_flag=123'.

Most importantly, invalid command line parameters are rejected with a suitable error message.

.. code-block:: shell-session

    $ cmdline/examples/ae108-CmdLineExample --enable_greeting=123
    the argument ('123') for option '--enable_greeting' is invalid. Valid choices are 'on|off', 'yes|no', '1|0' and 'true|false'

Using the Example
^^^^^^^^^^^^^^^^^

The full source code of the example is available in ``cmdline/examples/Example.cc``.
Here's what it looks like:

.. literalinclude:: cmdline-example.cc
    :language: cpp
    :linenos:

If you want to build it and try it out, then compile the executable target ``ae108-CmdLineExample`` and run it with command line options of your choice.

Outlook
-------

We've seen many of the features in action, but there's a bit more to explore.
For instance, it's possible to chain more than one call to ``withOption``:

.. code-block:: cpp

    cmdline::CommandLineOptionParser(std::cerr)
        .withOption(/* ... */)
        .withOption(/* ... */)
        .parse(argc, argv);

Also, there is another overload of ``withOption`` that permits to add flags without a help text.

You can find additional information in the API documentation of the library.

.. toctree::

    cmdline-src.rst

In addition, the tests in ``cmdline/test/CommandLineOptionParser_Test.cc`` showcase the features in common use cases.
