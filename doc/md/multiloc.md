# Implementing multilocus simulations

In order to simulate discontiguous genomic segments, one has two options within __fwdpp__:

1. Write mutation, recombination, and fitness policies that "do the right thing" for your model.  For example, a mutation policy would need to know about the mutation rate at each locus, and appropriately assign mutations with the correct positions, fitness effects, etc.
2. Write separate mutation and recombination policies for each locus, and a fitness policy that calculates the fitness of a diploid over all loci.

This document covers the latter method, which I call the "multilocus" part of __fwdpp__.    I won't give any examples of the former method, as I'm opposed to the idea of having to develop, debug, and mainting large complex policies.  But if you want to do things that way, read the [tutorial on policies](@ref md_md_policies), because what you are looking for is all in there.

# Simple policies using the multi-locus machinery

__fwdpp__ only supports the implementation of individual-based multilocus simulations, and it is unlikely that these features will be added to the gamete-based part of the library.

The main conceptual difference between this part of the library and the examples shown in the [tutorial on policies](@ref md_md_policies) is the following:

* Instead of a single mutation model policy, you implement one mutation model per "locus".  These policies are stored in a vector and passed to KTfwd::sample_diploid.
* Similarly, you implement a recombination policy per locus, and pass a vector of those policies along to KTfwd::sample_diploid.
* Instead of a single list of gametes, you have a vector of lists of gametes, where each list represents the current gametes at a particular locus.
* A diploid is now represented as a vector of pairs of iterators derived from the vector of lists of gametes.
* A fitness policy calculates individual fitnesses from that vector of pairs of iterators.

At this point, it may be most useful to look at a concrete example.  The program diploid_ind_2locus.cc is distributed with the library source code, and we'll break down its essential parts in the next few sections.

### A multilocus mutation model

### Separate recombination policies per locus

### A (trivial) mutilocus fitness model

This needs special mention

## Mechanics of the multilocus recombination

Document the logic here
