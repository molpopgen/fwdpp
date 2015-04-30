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

### Notes

* Remember to _properly_ initialize the base class from your derived mutation classes!
* All fitness models provided by the library require a mutation type to contain a double called \f$s\f$ representing the selection coefficient (or effect size).  For example, see KTfwd::site_dependent_fitness.  Hoewever, if you write your own fitness models, then your mutations can contain whatever they want. 

## Mutation lists

Mutations must be stored in double-linked lists.  In the C++ standard library, use std::list:

~~~{.cpp}
#include <list>
#include <fwdpp/diploid.hh>

std::list< KTfwd::mutation > mlist;
~~~

The standard list type works well, and simulations will be pretty fast.  However, we can do better.  In any forward simulations, mutations enter the population and then quickly exit because they are lost due to drift and/or selection.  Situations where small objects are coming and going rapidly screams for a memory pool allocation model, rather than the standard allocator that the above code block is implying:

~~~{.cpp}
#include <list>
#include <fwdpp/diploid.hh>

//This is equivalent to the above block,
//which relies on the default allocator
//being the standard one.
std::list< KTfwd::mutation, std::allocator<KTfwd::mutation> > mlist;
~~~

The [boost](http://www.boost.org) C++ libraries provide nicely-implemented memory pool objects:

~~~{.cpp}
#include <boost/container/list.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <fwdpp/diploid.hh>

//We see that typedefs lead to improved sanity:
using mutation_t = KTfwd::mutation;
using mutation_t_allocator = boost::pool_allocator<mutation_t>;
boost::container::list< mutation_t, mutation_t_allocator > mlist;
~~~

For simulations with large mutation and/or recombination rates, using the boost types can result is substantially faster simulations.  It is well worth it to have boost available on your systems!

Note that you can use the boost allocators with the standard containers, etc.

Further note that I have done _zero_ experimentation with parameters to the pool_allocator template, which itself can be tweaked in various ways.  

## Gamete types

The basic gamete type is KTfwd::gamete_base, which is a template type.  Fundamentally, gametes contain vectors if _iterators_ derived from mutation lists.  There is one such vector for "neutral" variants, another for "selected" variants.  

In order to instantiate a gamete template, you must provide the mutation type at a minimum:

~~~{.cpp}
#include <list>
#include <fwdpp/diploid.hh>

using mutation_t = KTfwd::mutation;
std::list< mutation_t > mlist;
//This defaults to assuming that std::list<mutation_t> is the mutation list type
using gamete_t = KTfwd::gamete_base<mutation_t>;
~~~

If your mutation list type is not std::list, then you must provide an additional template parameter:

~~~{.cpp}
#include <boost/container/list.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <fwdpp/diploid.hh>

//We see that typedefs lead to improved sanity:
using mutation_t = KTfwd::mutation;
using mutation_t_allocator = boost::pool_allocator<mutation_t>;
using mlist_t = boost::container::list< mutation_t, mutation_t_allocator >;
mlist_t mlist;
using gamete_t = KTfwd::gamete_base<mutation_t,mlist_t>;
~~~

That's basically it for gametes. 

### Notes

* The interface to Ktfwd::gamete_base is a little funny.  Ideally, you should be able to provide only the mutation list type, _e.g._,

~~~{.cpp}
using gamete_t = KTfwd::gamete_base<mutation_t,mlist_t>;
~~~

The gamete_t should be able to deduce the mutation_t from the mlist_t::value_type.  That's all true, but I first wrote the gamete type based around the idea of the mutation and the default std::list<mutation_t>.  Changing it now would break all existing simulations using the library.

* As of __fwdpp__ 0.3.0, KTfwd::gamete_base takes a third template argument, which currently defaults to KTfwd::tags::standard_gamete.  This tag has no effect on anything, and the library remains source-compatible with previous versions.  This "tag" is a placeholder for potential future library features making use of tag-dispatch methods for implementing more complex models.

## Gamete lists

Like mutations, gametes are stored in doubly-linked lists.  All of the stuff from the section on mutation lists applies here.

## Diploids

In fwdpp, a diploid is a pair of iterator derived from a gamete list, _e.g._,

~~~{.cpp}
using glist_t = std::list<gamete_t>;
using diploid_t = std::pair< glist_t::iterator, glist_t::iterator >;
~~~

Then, the entire population is a vector of diploids:

~~~{.cpp}
std::vector< diploid_t > diploids;
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
