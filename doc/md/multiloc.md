# Tutorial 2: Implementing multilocus simulations

In order to simulate discontiguous genomic segments, one has two options within __fwdpp__:

1. Write mutation, recombination, and fitness policies that "do the right thing" for your model.  For example, a mutation policy would need to know about the mutation rate at each locus, and appropriately assign mutations with the correct positions, fitness effects, etc.
2. Write separate mutation and recombination policies for each locus, and a fitness policy that calculates the fitness of a diploid over all loci.

This document covers the latter method, which I call the "multilocus" part of __fwdpp__.    I won't give any examples of the former method, as I'm opposed to the idea of having to develop, debug, and mainting large complex policies.  But if you want to do things that way, read the [tutorial on policies](@ref md_md_policies), because what you are looking for is all in there.

# Differences from a "single-locus" simulation.

The multi-locus API allows the simulation of non-contiguous, partially-linked regions.  The recombination between these
regions may be modeled as a Poisson process (_e.g._, the number of crossovers between regions $i$ and $i+1$ is a Poisson
random variable), or a binomial process (_e.g._, representing the traditional notion of a genetic distance).

In terms of data types, a diploid is now a vector of objects that, at minimum, provide the same public interface as
std::pair<std::size_t,std:size_t>.  

Mutation models are per locus/region, and are stored in a vector-like container.  The same is true of within- and
between-region recombination models.  For a system of $k$ loci, there must be $k$ mutations models, $k$ within-locus
recombination models, and $k-1$ between-locus recombination rates.

For concrete examples, see the following files:

* diploid_ind_2locus.cc, which is an example program
* sugar_multilocusTest.cc, which is part of the library's testing suite

# Limitations

The overload of fwdpp::sample_diploid for multi-region simulations is limited to modeling inter-locus recombination as
either Poisson or binomial, and not as a mixture.  This will addressed in a future version of the library.
