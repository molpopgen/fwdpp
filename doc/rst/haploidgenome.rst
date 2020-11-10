.. _haploidgenomes:

Haploid genomes
=====================================================

This section builds on the information found in :ref:`mutations`.

A haploid genome is the information inherited from a single parent.
For example, a diploid individual contains two haploid genomes.
In the case of obligate outcrossers, one genome may be labelled "maternal" and the other "paternal".

The type :type:`fwdpp::haploid_genome` represents a haploid genome in ``fwdpp``.
This type is itself a ``typedef`` of `fwdpp::haploid_genome_base`.

The key data members are:

* :member:`fwdpp::haploid_genome_base::n` is the number of times this genome is present in the simulation.
* :member:`fwdpp::haploid_genome_base::mutations` contains the indexes of mutations that do **not** affect genetic values and/or fitness.
* :member:`fwdpp::haploid_genome_base::smutations` contains the indexes of mutations that **do** affect genetic values and/or fitness.

The two classes of mutation keys are stored in different containers so that we can skip "neutral" mutations when calculating genetic values.

The type of the index containers is :type:`fwdpp::haploid_genome_base::mutation_container`.
The ``value_type`` of this container is :type:`fwdpp::uint_t`.

.. note::

    The field :member:`fwdpp::haploid_genome_base::n` is not equivalent to the frequency of a genome.
    Rather, it is the number of occurrences of a *specific* genome object.

Containers of genomes
+++++++++++++++++++++++++++++++++++++++++++++++++++++

In a simulation, genome objects must be stored in random-access containers of objects.

In other words, the ``value_type`` of your container must equal a genome type;

For example:

.. code-block:: cpp

   using genome_container = std::vector<fwdpp::haploid_genome>;
