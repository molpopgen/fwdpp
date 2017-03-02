# Overview of fwdpp features

This document is a quick tour of __fwdpp__'s features.  It also points out what is missing from the library that may or may not be addressed in the future.

The document assumes a working knowledge of C++11 and skips details, such as class constructor implementation, that are C++11 issues rather than __fwdpp__ issues.

Sub-sections named "aside" may be skipped--they contain some technical details that hopefully address questions that proficient C++ programmers may have.

## Library headers

The library is header-only, meaning that there is no runtime library to link to. Rather, the vast majority of library code is implemented as generic templates.

There are two main headers to be aware of:

~~~{.cpp}
//This includes the low-level library:
#include <fwdpp/diploid.hh>
//This includes the higher-level objects
//and functions:
#include <fwdpp/sugar/sugar.hpp>
~~~

## Library organization
 
The two headers correspond to two related parts of __fwdpp__.  The header `fwdpp/diploid.hh` contains low-level functionality.  This code base represents how the library grew over time.  

While it is possible to write your simulations using only the part of __fwdpp__, doing so is more cumbersome than it needs to be. During our work in developing simulations for our research, a series of best practices have evolved, which developed into the "sugar" part of the library (`fwdpp/sugar/sugar.hpp`).  

We recommend using "sugar" features whenever possible.  They call the low-level functions an often have simpler APIs.  In the sections below, function calls from "sugar" will be noted.

Both the low-level and sugar parts of __fwdpp__ are in the namespace `KTfwd`.

## Built-in types

This section describes the low-level types required by a simulation.

### An unsigned integer

The library defines this type:

~~~{.cpp}
#include <cstdint>
namespace KTfwd
{
    using uint_t = std::uint32_t;
}
~~~

Thus, the name `KTfwd::uint_t` refers to a 32-bit, unsigned integer.

### Diploids

When modeling a single contiguous genomic segment, a diploid is simply a pair of gametes.  The way to represent that using __fwdpp__ is actually with a type from the C++ standard library:

```cpp
#include <utility>
#include <cstdint>

// This is a typedef for a diploid
using diploid = std::pair<std::size_t,std::size_t>
```

#### Custom diploids

Sometimes, you need more information associated with your diploid type.  You can define a custom diploid provided that they provide the same _minimal API_ the pair from the previous section.

~~~{.cpp}
/* This is a custom diploid
 * that is essentially a hand-coded
 * instantation of pair<size_t,size_t>
 */
struct minimal_custom_diploid
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first,second;

    //define constructors, etc., as needed
};
~~~

The code block above will be considered a custom diploid type by __fwdpp__.  Of course, this one is trivial.  Realistic examples may include an integer type representing biological sex, for example:

~~~{.cpp}
struct diploid_with_sex
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first,second;
    int sex;

    //define appropriate constructors, etc.
};
~~~

#### An aside

Why does __fwdpp__ not simply have custom diploid types inherit from the standard pair type?  Well, the standard library pair does not have a virtual destructor, meaning that it cannot be used as a base type.  

#### Possible future changes

Future versions of the library may add a custom diploid base type.  However, diploids get copied/moved with some regularity, meaning that it is helpful for them to be trivally constructible, etc., which discourages a class hierarchy.

### Multi-locus/multi-region diploids

When modeling multiple contiguous regions separated by some genetic distance, a diploid is represented as a vector of single-region diploid types.  For example:

~~~{.cpp}
#include <vector>

using multilocus_diploid = std::vector<diploid>;
~~~

You can also have vectors of custom diploids:

~~~{.cpp}
#include <vector>

using multilocus_custom_diploid = std::vector<diploid_with_sex>;
~~~

#### Possible future changes

When using a custom diploid, the `vector` will, of course, copy any of the associated information, leading to redundancy and a bit more memory consumption.  In practice, we only update the information for the first element of the vector. However, a future version of the library may define a new data type to prevent the redundancy.

### Gametes

A gamete is defined in the low-level part of the library (`fwdpp/diploid.hh`)`.

A gamete is quite simple in structure, containing the following:

* An unsigned integer representing how many times this specific instance of the gamete exists in the population.
* A container of "keys" to neutral mutations
* A container of "keys" to selected mutations

A gamete is quite simple in structur, containing the following:

* An unsigned integer representing how many times this specific instance of the gamete exists in the population.
* A container of "keys" to neutral mutations
* A container of "keys" to selected mutations
* The key type is `std::uint32_t`.

Let's go right to its declaration:

~~~{.cpp}
namespace KTfwd
{
    template <typename TAG = tags::standard_gamete>
    struct gamete_base
    {
        //Count in population
        uint_t n;
        //Dispatch tag type
        using gamete_tag = TAG;
        using index_t = std::uint32_t;
        using mutation_container = std::vector<index_t>;
        //Container of mutations not affecting
        //trait value/fitness ("neutral mutations")
        mutation_container mutations;
        //Container of mutations affecting
        trait value/fitness ("selected mutations")
        mutation_container smutations;

        //Constructor
        gamete_base(const uint_t &icount) noexcept;

        //Perfect-forwarding constructor
        template <typename T>
        gamete_base(const uint_t &icount, T &&n, T &&s) noexcept;
        
        //Destructor is virtual,
        //so you may inherit from this type
        virtual ~gamete_base() noexcept {}

        //Copy/move constructors, etc.
        //are omitted for space. In C++11
        //lingo, they are declared "default"

        //Equality comparison
        inline bool
        operator==(const gamete_base<TAG> &rhs) const;
    };
}
~~~

### Mutations

A mutation is an object that publicly inherits from `KTfwd::mutation_base`, which as the following declaration:

~~~{.cpp}
namespace KTfwd
{
    struct mutation_base
    {
        //Mutation position
        double pos;
        /*
          16 bits of extra data to be associated w/this type.
          Do with it what you will. Fits into padded space in this struct,
          and doesn't affect sizeof(mutation).
        */
        std::uint16_t xtra;
        //Is the mutation neutral or not?
        bool neutral;
        //Constructor
        //(Be sure to call this from 
        //derived classes!)
        mutation_base(const double &position,
                      const bool &isneutral = true,
                      const std::uint16_t x = 0) noexcept;
        //destructor
        virtual ~mutation_base() noexcept {}
        //The usual suite of constructors
        //are declared "default", and omitted
        //here.
};
}
~~~

The base mutation class is very simple.  It records a position and the "neutrality" of the variant.  In fact, it is too simple to be useful.  The "sugar" layer of the library provides a set of more generally useful mutation types.  In practice, we find that these satisfy most needs:

* `KTfwd::popgenmut` is a mutation with an effect size (`s`), dominance (`h`) and a record of when it appeared (`g`).  Note that __fwdpp__ does not care if `s` is an effect size on a phenotype or a "selection coefficient" in the way that population geneticists typically define the term.
* `KTfwd::generalmut` where `s` and `h` are stored as `std::array<,doublestd::size_t>`, allowing for multiple `s/h` values to be associated with a variant.  Think of using this array to simulate pleiotropic effect sizes, for example.
* `KTfwd::generalmut_vec` is largely equivalent to `KTfwd::generalmut`, but `std::vector<double>` replaces the `std::array`.  Possible use cases could involve different variants having different numbers of pleitropic effects.  Or, and perhaps more practically, this type can be used in other language environments that do not understand `std::array` (think [Cython](http://www.cython.org)).

The three types listed above are included via `#include <fwdpp/sugar/sugar.hpp>`.  See the [reference manual](http://molpopgen.github.io/fwdpp/doc/html/index.html) for more details about them.

### Mutations + Gametes + Diploids = (almost) a population.

Let's put the previous three sections into context:

* We will need a container of mutations for our simulation.  For example, `using mcont_t = std::vector<KTfwd::popgenmut>;`cpp.

#### Aside: `std::uint32_t` for mutation keys?

The astute reader will notice that `KTfwd::gamete_base` stores mutation keys as `std::uint32_t` while `std::vector<any mutation type>::size_type` is almost certainly `std::size_t`, which will be 64 bits on most systems.  In practice, `sizeof(a mutation type)` is sufficiently large that it is very unlikely that anyone can store more than $2^32-1$ such objects.  Thus, __fwdpp__ intentionally makes a compromise on the integer width.  A nice side-effect is that we save a lot of memory!  Those keys end up being a lot of the storage in a simulation of a large population/large genomic region.

## Representing a population

## Modeling the biology

### Mutation

### Recombination

### Fitness
