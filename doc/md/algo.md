# Algorithmic details of recombination

(This document first appeared in __fwdpp__ 0.3.3.)

The release of __fwdpp__ 0.3.3 contained two significant improvements in the scalability of the code for large simulations.

## Simplifying crossing-over

Through __fwdpp__ 0.3.2, the algorithm for crosssing over went like this:

* Take the two gametes in a parent
* Recombine them at a set of breakpoints determined by the user-defined policy
* Randomly-pick one of the two recombinant gametes to pass on to the offspring.

In "pictures", the algorithm went like:

* Starting with parental gametes: AAAAAA and BBBBBB
* Ending with recombinant gametes: AAABBB and BBBAAAA
* Pick one to pass on (say, AAABBB), and discard the other)

On obvious problem with this approach is that we are creating a gamete that we don't do anything with: the one that doesn't get passed on simply gets discarded, and is thus a waste of all the copying (of pointers to mutations) that it took to create that gamete.

Version 0.3.3 of the library introduced the following change (which I honestly wish I'd thought of several years ago...):

* Take the two gametes in a parent, pg1 and pg2
* Randomly pick one to be the descendant.  If it is pg2, swap the pointers to pg1 and pg2.
* Create a single recombinant and assign the appropriate pointer as the new value of pg1

In "pictures":

* Starting with parental gametes: pg1 = AAAAAA and pg2 = BBBBBB
* Half the time, swap pg1 and pg2
* Generate a single recombinant, either AAABBB _or_ BBBAAA

What are the consequences of this new method?  Let's start with the good:

* Reduced run-time.  I _think_ the benefit is typically 15-20%, but this change was not "good science", as I also changed how memory for storing recombinant gametes was allocated (see next subsection).

And now, the bad:

* For a given random number seed, the output is different for programs compiled using 0.3.2 and 0.3.3.  This is because the "Mendel" step occurs earlier in 0.3.3, thus changing the context of all subsequent calls to the random number generator.

### Changing memory allocation patterns

Through __fwdpp__ 0.3.2, the library would allocate memory for the two recombinant gametes.  This occurred for every recombination event, meaning that for \f$\rho = 4Nr\f$, the expected number of containers created each generation was \f$4\rho\f$, as there are \f$4Nr\f$ expected recombination events each generation (2 parents per 2N offspring, each parent recombining at \f$r\f$ positions on average), and 4 containers needed per recombining parent (because pointers to neutral and selected mutations are stored separately).

The method describe above worked ok, but ultimately results in a massive number of calls for new memory over the course of a simulation.

In 0.3.3, we can cut the number of allocations in half because we generate half as many recombinant chromosomes (see above).   However, the library takes the following steps to reduce allocations even further:

* Programs must now manage the allocation of the two vectors required to store a recombinant.  Let's call them "neutral" and "selected".
* They are declared in main(), and the programmer should use "reserve" to pre-allocate some memory for them.
* These containers are passed into recombination policies and emptied/filled as required.

The positives of the new approach are:

* Reduced  memory use.  We can now just let the vector class adjust its capacity as needed, and grow as needed during the simulation, resulting in many fewer requests for memory.
* In combo with the new crossing-over method (see above), run times are reduced.

The negative:

* Passing these containers to recombination policies represents an API change.

## Checking that a recombinant gamete is unique

__fwdpp__ uses doubly-linked lists for storing pointers (iterators) to mutations and gametes.  The advantage of these lists is that pointers remain valid after insertion/deletion to/from lists. (This is in contrast to the case with a vector/array, where insertion/removal triggers a reallocation, which often moves memory, resulting in "pointer/iterator invalidation".)   The down-side of lists is that moving through them is relatively expensive.  They fragment memory over time because they use discontinuous storage, and therefore searching through them is expensive.

One of the goals of the library is to store objects once and only once.  Thus, after a recombination, we need to do the following:

* Ask if the recombinant gamete is unique.
* If not, we simply assign the pointer to the recombinant as the pointer to the version that we already have in a list
* If so, we insert the new gamete at the end of the list, and store a pointer to it.

Through __fwdpp__ 0.3.2, the library searched the entire gamete list to see if the recombinant already exists.   This search took \f$O(g)\f$ operations, where \f$g\f$ is the number of unique gametes in the population.  Such linear searches are slow for linked lists and scale terribly with increasing population size and/or mutation and/or recombination rate.

In 0.3.3, the linear searches are replaced with \f$O(log(g))\f$ lookups via a lookup table that is calculated internally every generation.

The lookup table associates the total number of mutations in a gamete (neutral + selected mutations) and the iterator to that gamete in the linked list.  Given a new gamete, we can find the _range_ of gametes in the lookup table with the same total number of mutations in logarithmic time.  If the size of that range is zero, then the new gamete is by definition unique, and we insert it.  Otherwise, we must compare the recombinant to all gamtes in the range, which is a linear serach, but this time over a much smaller range than what __fwdpp__ 0.3.2 was doing.

The plus sides:

* Massive improvements in the scalability of the simulation.  For the cost of building the lookup table, we replace \f$O(\rho/2)\f$ linear searches of a linked list with \f$O(\rho/2)\f$ logarithmic-time lookups.  Further, _we rarely have to do the secondary linear search!_

The down side (quibble):

* We must pass the lookup table on to the recombination policy, and thus we need an API change.   However, we've already changed the API in 0.3.3 because of the changes to how recombination works, so this is a very nit-picking down-side.

### Details on the lookup table

(This section exists in hope that KRT reads it before trying to reinvent the wheel in the future.  Fat chance of that...)

During testing, I tried the following schemes for the lookup tables:

* std::vector< std::pair<std::int32_t, glist::iterator> >, sorted using either std::stable_sort or std::sort.
* std::multimap< std::int32_t, glist::iterator >
* std::unordered_multimap< std::int32_t, glist::iterator >

I benchmarked each approach with the following command line: 

~~~{sh}
#Keep seed the same for all three...
diploid_ind 10000 4000 4000 100000 10 1 SEED
~~~

All three lookup tables had very similar run-times and peak RAM use.  However, (my implementation using) the unordered_multimap resulted in simulations with incorrect distributions of summary statistics.

The library currently defaults to using the std::multimap.  Programs may be compiled using "sorted vector of pairs" approach by passing -DFWDPP_VECTOR_GLOOKUP to the compiler/preprocessor.
