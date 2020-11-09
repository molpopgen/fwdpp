Compiling, installing, etc.
===========================================================

Obtaining the source
++++++++++++++++++++++++++++++++++++++++++++++++++

The source code can be obtained from ``GitHub``.
For example,

.. code-block:: bash

   git clone https://github.com/molpopgen/fwdpp
   cd fwdpp
   git submodule update --init --recursive

.. note::

    It is important to update the submodules!
    Some of the examples and several of the unit
    tests make use of ``tskit``, which is one
    of the submodules.

.. warning::

    Be careful if you download the source code from the "releases" section of ``GitHub``.
    Do not download the files that ``GitHub`` automatically generates!
    They do not include the submodule.

Dependencies
++++++++++++++++++++++++++++++++++++++++++++++++++

1. A compiler that supports ``C++14`` or later.
2. GNU autotools
3. The GNU Scientific Library
4. The boost libraries.
   The ``boost-test`` and ``boost-program-options`` are used for the test suite and for some examples programs, respectively.

All of these dependencies are available via package managers.
On Debian/Ubuntu-like systems, you may use `apt` to install them.
On other operating systems, ``conda`` and/or ``brew`` will do the job.

Building examples and running tests
++++++++++++++++++++++++++++++++++++++++++++++++++

.. note::

    The test suite will not be compiled if the boost test library is not found on your system.

To build the examples and the test suite:

.. code-blocK:: bash

   autoreconf --install
   ./configure
   make

The compilation can take a long time.
To enable parallel compilation, you may execute ``make -j 6`` to use six threads, for example.

To run the tests:

.. code-block:: bash

   make check

To install the library:

.. code-block:: bash

   # This may require sudo
   make install


Changing the compilation flags
---------------------------------------------

The default compilation flags are ``-DNDEBUG -O2 -g``.
The ``-DNDEBUG`` specifies compiling in "release mode".
To compile in "debug" mode:

.. code-block:: bash

   ./configure --enable-debug=yes

.. note::

    Compiling in debug mode makes simulations run a lot slower!
    The test suite is always compiled in debug mode.

To change where the library will be installed:

.. code-block::bash

   # This will install to ~/include with a "make install"
   ./configure --prefix=$HOME

To change the optimization levels, etc.:

.. code-block:: bash

   CXXFLAGS="-O3 -W -Wall -Wconversion" ./configure

The above example enables more aggressive optimizations and sets more warning flags.

The default C++ mode is ``C++14``.
To enable `C++17`:

.. code-block:: bash

   ./configure --enable-cpp17=yes

If you have dependencies installed in non-standard locations, then you may need to provide those locations to the configure script.
When doing so, keep in mind that ``CPPFLAGS`` is used to change where the compiler looks for headers, as it is the variable affecting the *preprocessor*.
The ``CXXFLAGS`` variable is for the *compiler* and **not** the preprocessor!
