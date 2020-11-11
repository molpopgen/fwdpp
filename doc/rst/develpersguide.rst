.. _developersguide:

Developers guide
================================================================

Calculating code coverage
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Code coverage requires `lcov`, which will be available via your favorite package manager.

The coverage is calculated from the test suite and is automated via a ``make`` target:

.. code-block:: bash

   make -C testsuite coverage-local

To clean up after a coverage calculation:

.. code-block:: bash

   make clean-local

The output will be in the directory ``fwdpp_coverage``, which will be generated in the root directory of the source repository.
To view the results, point your browser to ``fwdpp_coverage/index.html``.

.. note::

   Coverage calculations are a bit odd for header-only libraries like this one.
   The output only applies to library files that are included in the test suite.
   Thus, there may be files for which no coverage calculation is possible.
   We are working to ensure that all files do get covered.
