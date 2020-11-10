.. _mutations:

Mutations
=============================================

Most types of simulations must have mutations in order to do anything "interesting".
``fwdpp`` allows you to define your own mutation type.

A valid mutation type publicly inherits from :class:`fwdpp::mutation_base`.
This class defines a very minimal interface.
It has no concept of an "effect size" for a mutation, etc..
To add the relevant concepts for your model, create a derived class.

For example:

.. literalinclude:: ../../examples/custom_mutation_example.hpp

Containers of mutations
+++++++++++++++++++++++++++++++++++++++++++++++++++++

In a simulation, mutation objects must be stored in random-access containers of objects.
You must not use, for example, smart pointers to :class:`fwdpp::mutation_base`.

In other words, the ``value_type`` of your container must equal a type derived from the mutation base class.

For example:

.. code-block:: cpp

   using mutation_container = std::vector<mutation>;

By using such a container, each mutation is represented only once in a simulation.
To represent individual genotypes, we need a method of tracking which mutations are present in genomes.
See :ref:`here <haploidgenomes>` for details.
