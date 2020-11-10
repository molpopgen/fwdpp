.. _diploids:

Diploids
===================================

This section builds on the information in :ref:`haploidgenomes`.

The simplest representation of a diploid is a pair of integers referring to the simulation's genome container.
``fwdpp`` treats ``std::pair`` with integer types as the simplest form of diploid.
At compile time, the library ensures that ``first_type`` and ``second_type`` are both integer types and that both are the same type.

For example:

.. literalinclude:: ../../examples/diploid_example.cc

In the above example, the static assertions tell us that our ``std::pair`` are diploids.
They also tell us that they are not custom diploids..

Diploids are more than just pairs of genomes.
You may want to record other meta data along with the genome indexes.
To do so, you may define your own types that duck-type the ``std::pair`` interface, adding whatever features you would like.
We refer to such types as custom diploid types:

.. literalinclude:: ../../examples/custom_diploid_example.cc
   :lines: 7-30

Internally, ``fwdpp`` doesn't care much about custom diploids versus ``std::pair``.
The support for custom diploids exists to give modeling flexibility.
However, the library is flexible enough to let you take other paths as well:

* You could represent meta data as a separate array of structures.
* You could represent the meta data as a separate structure of arrays.

Containers of diploids
+++++++++++++++++++++++++++++++++++++++++++++++++++++

In a simulation, diploid objects must be stored in random-access containers of objects.

In other words, the ``value_type`` of your container must equal a diploid type;

For example:

.. code-block:: cpp

   using diploid_container = std::vector<my_diploid_type>;
