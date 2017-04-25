# Tutorial 0: Data types in fwdpp

This document introduces the fundamental data types for simulations implemented using __fwdpp__.

The relevant header file defining these types is 

~~~{.cpp}
#include <fwdpp/forward_types.hpp>
~~~

This header will be included when you include the main library header via

~~~{.cpp}
#include <fwdpp/diploid.hh>
~~~

__Note:__ this document only covers single-population, individual-based simulations in detail.

## Mutation types

The most basic type determining the behavior of a simulation is a mutation.  Any __fwdpp__-based simulation must define a mutation type that publicy inherits from KTfwd::mutation_base.  This base class contains several public data types:

~~~{.cpp}
//The mutation position
mutable double pos;
//The number of occurrences in the population
unsigned n;
//Does this mutation affect fitness/trait value?
bool neutral;
//This is used and updated internally by other functions w/in the library
bool checked;
~~~

This type is usable in a simulation, but not in an interesting simulation, as there is no selection coefficient, etc.  In order to add things to a mutation type, you extend KTfwd::mutation_base via public inheritance.  You may find an example of that in the next tutorial (@ref md_md_policies) and in the following types provided by the library:

* KTfwd::mutation
* KTfwd::popgenmut
* KTfwd::generalmut
* KTfwd::generalmut_vec

### Notes

* Remember to _properly_ initialize the base class from your derived mutation classes!
* All fitness models provided by the library require a mutation type to contain a double called \f$s\f$ representing the selection coefficient (or effect size).  For example, see KTfwd::site_dependent_genetic_value.  Hoewever, if you write your own fitness models, then your mutations can contain whatever they want. 

## Mutation containers

Mutations must be stored in a container supporting random access via operator[].  Typically, this will be std::vector.
For example: 

~~~{.cpp}
#include <vector>
#include <fwdpp/diploid.hh>

//Create a typedef for a container of mutations:
using mcont_t = std::vector<KTfwd::mutation>;
~~~

## Gamete types

The basic gamete type is KTfwd::gamete_base, which is a template type. A gamete contains two std::vector<std::size_t>
called _mutations_ and _smutations_.  The values contained in these vectors correspond to the mutations present in the
gamete, and are used as indexes into the mutation container described above.  The containers _mutations_ and
_smutations_ contain indexes (or "keys") to neutral and non-neutral mutations, respectively.

A gamete also contains an unsigned integer, _n_, representing how many times the gamete is present in the simulation.

__NOTE:__ __fwdpp__ makes no attempt to collect identical gametes together. In other words, imagine a recombination
event results in a gamete identical to one that already exists.  The library will _not_ check for this, because the
comparison is relatively slow.  Thus, the _n_ value will tend towards 1 over time, particularly for simulations with
large mutation/recombination rates.

That's basically it for gametes. 


## Gamete containers

Like mutations, gametes are stored in vector-like containers.   All of the stuff from the section on mutation containers applies here.

## Diploids

In fwdpp, the simplest diploid is std::pair<std::size_t,std::size_t>.  These two integers are indexes/keys referring to
the diploid's two gametes.

The information represented in a diploid may be augmented via custom diploid types. See @ref md_md_customdip for
details.

A simulation of a single deme makes use of a vector of diploids:

~~~{.cpp}
std::vector<std:pair<std::size_t,std::size_t> > diploids;
~~~

A simulation of multiple demes makes use of a vector of vector of diploids:

~~~{.cpp}
using diploid_t = std::pair<std::size_t,std::size_t>;
using dipvector_t = std::vector<diploid_t>;
std::vector<dipvector_t> diploids;
~~~

## Accessing the data in a population

With the above definitions, we are ready to access all of the information in a population.  Here's a quick run-through:

~~~{.cpp}
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/popgenmut.hpp>
#include <iostream>

using mtype = KTfwd::popgenmut;
using mcont_t = std::vector<mtype>;
using gamete_t = KTfwd::gamete;
using gcont_t = std::vector<gamete_t>;
using diploid_t = std::pair<std::size_t,std::size_t>;
using dipcont_t = std::vector<diploid_t>;

mcont_t mutations;
gcont_t gametes;
dipcont_t diploids;

//Run a simulation, etc., now...

//Go over all diploids, printing out the neutral and selected mutation
//positions of each variant in each diploid.
void print_mutations(const gamete_t::mutation_container & mc, const mcont_t & mutations)
{
    for(const auto & key : mc) {
    std::cout << mutations[key].pos << ' ';
    }
    std::cout << '\n';
}

for(const auto & dip : diploids)
{
   //neutral mutations, first gamete
   print_mutations(gametes[dip.first].mutations,mutations); 
   //selected mutations, first gamete
   print_mutations(gametes[dip.first].smutations,mutations);
   //Repeat process for second gamete
   print_mutations(gametes[dip.second].mutations,mutations);   
   print_mutations(gametes[dip.second].smutations,mutations);   
}

~~~


## What next

A simulation is implemented as a set of _policies_ that interact with the above data types:

* @ref md_md_policies

__fwdpp__ 0.3.0 introduced a "sugar" layer of code intended to make simulation development easier:

* @ref md_md_sugar

The examples are the best place to look for how simulations are actually implemented:

* @ref md_md_examples

There is some (limited) support for quick-starting projects based on __fwdpp__:

* @ref md_md_devtools

Finally, there are more advanced topics:

* @ref md_md_multiloc
* @ref md_md_serialization
