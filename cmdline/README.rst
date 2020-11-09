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

.. code-block:: cpp

    #include <ae108/cmdline/CommandLineOptionParser.h>
    #include <iostream>

    int main(const int argc, const char *const *const argv) {
        /**
         * @todo use command line options
         */
    }

Construction
^^^^^^^^^^^^

To parse the command line parameters, we'll construct a ``CommandLineOptionParser`` by providing a stream.
Let's say we want to print error messages and warnings to ``stderr``.

.. code-block:: cpp

    namespace cmdline = ae108::cmdline;
    cmdline::CommandLineOptionParser(std::cerr)
        /**
         * @todo define available options
         */
        ;

Adding Options
^^^^^^^^^^^^^^

The only thing that is missing are the options that we would like to parse.
In this example, we are going greet the world if greetings are enabled via a command line flag.
The ``CommandLineOptionParser`` class provides a ``withOption`` method that we are going to use to achieve that.

In this example we'll add a flag ``--enable_greeting`` (with a short form ``-g``) together with a help text describing the flag.

.. code-block:: cpp

    auto enable_greeting = bool{false};
    cmdline::CommandLineOptionParser(std::cerr)
        .withOption("enable_greeting,g", "Print a greeting.", &enable_greeting)
        .parse(argc, argv);

    if (enable_greeting) {
        std::cout << "Hello world!" << std::endl;
    }

Let's see it in action.

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

Example Source Code
^^^^^^^^^^^^^^^^^^^

The full source code of the example is available in ``cmdline/examples/Example.cc``.
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
