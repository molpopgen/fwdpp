# Overview of fwdpp features

This document is a quick tour of __fwdpp__'s features.  It also points out what is missing from the library that may or may not be addressed in the future.

The document assumes a working knowledge of C++11 and skips details, such as class constructor implementation, that are C++11 issues rather than __fwdpp__ issues.

Sub-sections named "aside" may be skipped--they contain some technical details that hopefully address questions that proficient C++ programmers may have.

## Library headers

The library is header-only, meaning that there is no runtime library to link to. Rather, the vast majority of library code is implemented as generic templates.

There are two main headers to be aware of:

```cpp
//This includes the low-level library:
#include <fwdpp/diploid.hh>
//This includes the higher-level objects
//and functions:
#include <fwdpp/sugar/sugar.hpp>
```

## Library organization
 
The two headers correspond to two related parts of __fwdpp__.  The header `fwdpp/diploid.hh` contains low-level functionality.  This code base represents how the library grew over time.  

While it is possible to write your simulations using only the part of __fwdpp__, doing so is more cumbersome than it needs to be. During our work in developing simulations for our research, a series of best practices have evolved, which developed into the "sugar" part of the library (`fwdpp/sugar/sugar.hpp`).  

We recommend using "sugar" features whenever possible.  They call the low-level functions an often have simpler APIs.  In the sections below, function calls from "sugar" will be noted.

Both the low-level and sugar parts of __fwdpp__ are in the namespace `KTfwd`.

## Built-in types

This section describes the low-level types required by a simulation.

### An unsigned integer

The library defines this type:

```cpp
#include <cstdint>
namespace KTfwd
{
    using uint_t = std::uint32_t;
}
```

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

```cpp
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
```

The code block above will be considered a custom diploid type by __fwdpp__.  Of course, this one is trivial.  Realistic examples may include an integer type representing biological sex, for example:

```cpp
struct diploid_with_sex
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first,second;
    int sex;

    //define appropriate constructors, etc.
};
```

#### An aside

Why does __fwdpp__ not simply have custom diploid types inherit from the standard pair type?  Well, the standard library pair does not have a virtual destructor, meaning that it cannot be used as a base type.  

#### Possible future changes

Future versions of the library may add a custom diploid base type.  However, diploids get copied/moved with some regularity, meaning that it is helpful for them to be trivally constructible, etc., which discourages a class hierarchy.

### Multi-locus/multi-region diploids

When modeling multiple contiguous regions separated by some genetic distance, a diploid is represented as a vector of single-region diploid types.  For example:

```cpp
#include <vector>

using multilocus_diploid = std::vector<diploid>;
```

You can also have vectors of custom diploids:

```cpp
#include <vector>

using multilocus_custom_diploid = std::vector<diploid_with_sex>;
```

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

```cpp
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
        //genetic value/fitness ("neutral mutations")
        mutation_container mutations;
        //Container of mutations affecting
        //genetic value/fitness ("selected mutations")
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
```

The library provides the following typedef for convenience:

```cpp
namespace KTfwd
{
    using gamete = gamete_base;
}
```

#### Aside: gamete_base's template type?

The "TAG" type, which defaults to `KTfwd::tags::standard_gamete` can be used for compile-time dispatch.  Currently, it is not used for anything, and it is only there in case of some (currently unknown) future need to distinguish different types of gametes at compile time.  The tag type is given a default type value so that it can be effectively ignored for right now.

### Mutations

A mutation is an object that publicly inherits from `KTfwd::mutation_base`, which as the following declaration:

```cpp
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
```

The base mutation class is very simple.  It records a position and the "neutrality" of the variant.  In fact, it is too simple to be useful.  The "sugar" layer of the library provides a set of more generally useful mutation types.  In practice, we find that these satisfy most needs:

* `KTfwd::popgenmut` is a mutation with an effect size (`s`), dominance (`h`) and a record of when it appeared (`g`).  Note that __fwdpp__ does not care if `s` is an effect size on a phenotype or a "selection coefficient" in the way that population geneticists typically define the term.
* `KTfwd::generalmut` where `s` and `h` are stored as `std::array<double,std::size_t>`, allowing for multiple `s/h` values to be associated with a variant.  Think of using this array to simulate pleiotropic effect sizes, for example. Or use a `std::array<double,2>` to have different effects in females vs males.  Etc.
* `KTfwd::generalmut_vec` is largely equivalent to `KTfwd::generalmut`, but `std::vector<double>` replaces the `std::array`.  Possible use cases could involve different variants having different numbers of pleitropic effects.  Or, and perhaps more practically, this type can be used in other language environments that do not understand `std::array` (think [Cython](http://www.cython.org), which does not currently support non-type template parameters required to use `std::array`).

The three types listed above are included via `#include <fwdpp/sugar/sugar.hpp>`.  See the [reference manual](http://molpopgen.github.io/fwdpp/doc/html/index.html) for more details about them.

If you have a need for a feature in a mutation that cannot be accomodated by one of the above types, simply define your own.  Your mutation type must publicly interit from `KTfwd::mutation_base`.

### Mutations + Gametes + Diploids = (almost) a population.

Let's put the previous three sections into context:

* We will need a container of mutations for our simulation.
* We will need a container of gametes for our simulation.
* We will need a container of diploids for our simulation.

Here, a container refers to a contiguous-memory container that is API compatible with the C++11 `std::vector`. In practice, we use the standard library types for all of our simulations.  Thus, we may define the following:

```cpp
#include <vector>
#include <utility> //for std::pair
#include <cstdlib> //for std::size_t
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar.hpp> //for KTfwd::popgenmut

using mcont_t = std::vector<KTfwd::popgenmut>;
using gcont_t = std::vector<KTfwd::gamete>;
using diploid = std::pair<std::size_t,std::size_t>;
using dipcont_t = std::vector<diploid>;
```

Let's think through these typedefs in light of the definitions of mutation, gamete, and diploid:

* A mutation is a fundamental type, and they are stored in a vector.
* Gametes are stored in vectors.  Gametes themselves contain vectors of `std::uint32_t`, which I said were "keys" to mutations.  Those keys represent the _indexes_ of the mutations present in a gamete.  "Neutral" and "selected" mutations are separated in order to speed up genetic value calculations.  Further, the keys are stored in _ascending_ order according to _mutation position_ (`KTfwd::mutation_base::pos`).  Storing things this way means that we can use _binary searches_ in algorithms like recombination that depend on mutation positions.
* A diploid is a pair of gametes.  The minimal diploid is `std::pair<std::size_t,std::size_t>`.  Those `size_t` objects are the _indexes_ into the gamete container.

Thus, a simulation requires three containers that all need to talk to each other:

```cpp
//store our mutations
mcont_t mutations;
//Store our gametes,
//which themselves contain
//location info referring
//to mutations
gcont_t gametes;
//Store diploids, which
//contain the location of 
//their gametes in gametes.
dipcont_t diploids;
```

The above typedefs and code block refer to a simulation of a single contiguous region and a single population/deme. For a simulation of multiple, partially-linked regions in a single deme, we can say:

```cpp
//rely on diploid def'n above
using multilocus_diploid = std::vector<diploid>;
using mlocus_dipcont_t = std::vector<multilocus_diploid>;
```

For a multi-population simulation of a single, contiguous genomic region, we can define:

```cpp
using metapop = std::vector<dipcont_t>;
```

If you're following along so far, you'll realize something odd.  The definitions `mlocus_dipcont_t` and `metapop` are the same type ! That is true--both are aliases for `std::vector<std::vector<diploid>>`. They are distinguished from one another in the library because operations on these different population types require different argument types.  In other words, template tricks are involved. 

#### Aside: `std::uint32_t` for mutation keys?

The astute reader will notice that `KTfwd::gamete_base` stores mutation keys as `std::uint32_t` while `std::vector<any mutation type>::size_type` is almost certainly `std::size_t`, which will be 64 bits on most systems.  In practice, `sizeof(a mutation type)` is sufficiently large that it is very unlikely that anyone can store more than $2^32-1$ such objects.  Thus, __fwdpp__ intentionally makes a compromise on the integer width.  A nice side-effect is that we save a lot of memory!  Those keys end up being a lot of the storage in a simulation of a large population/large genomic region.

## Representing a population

Sadly, we need more than containers of mutations, gametes and diploids.  There's more to book-keep, requiring additional data structures:

* We need a vector keeping track of mutation frequencies.
* It would be handy to have a vector of fixed mutations, so that we can remove such mutations from our population (for models where doing so would be appropriate).
* If we are tracking fixations, we probably want to track fixation times, leading to another vector recording those.
* When simulating infinitely-many sites mutation models, an efficient lookup table of current mutation positions would be a good idea.

This is where the "sugar" part of __fwdpp__ comes in handy.  While it would be possible to manually define all of these additional container types, it would be better if we could encapsulate all of these concepts into classes. The sugar layer defines three template types to help you: 

* `KTfwd::singlepop` to represent a contiguous region and a single deme
* `KTfwd::metapop` to represent a contiguous region and multiple demes
* `KTfwd::multiloc` to represent multiple partially-linked regions and a single deme.

Each of the above takes two type parameters as template arguments.  These are a mutation type and a diploid type, respectively.  Further, the diploid type defaults to `std::pair<std::size_t,std::size_t>`.  For example:

```cpp
#include <fwdpp/sugar.hpp>

//A population type where 
//the mutation type is
//KTfwd::popgenmut and
//a diploid is pair<size_t,size_t>
using singlepop_t 
    = KTfwd::singlepop<KTfwd::popgenmut>;

//A different population type based
//around a custom diploid
using singlepop_custom_t 
    = KTfwd::singlepop<KTfwd::popgenmut,diploid_with_sex>;
```

The rest of the types are filled in automatically for you.  These types always use `std::vector` for contiguous-memory containers.  The lookup type for mutation positions is defined in terms of `std::unordered_set`.

The `metapop` and `multiloc` template types work similarly.

These containers are defined as C++11 template aliases.  The `singlepop`, `metapop`, and `multiloc` types are template aliases for `KTfwd::sugar::singlepop`, `KTfwd::sugar::metapop`, and `KTfwd::sugar::multiloc`, respectively.  The types in namespace `KTfwd::sugar` take more type names as template parameters.  These types are the various container types, etc.  Thus, if you wish to use a vector type other than `std::vector` (perhaps `std::vector` with a custom allocator, or `boost::vector`), then you may define a new template alias in terms of those container types and it will "just work", provided that the container type's API matches those of the C++ standard library types. 

At this point, I'll refer you to the [reference manual](http://molpopgen.github.io/fwdpp/doc/html/index.html) for more detail on these types. 

## Recycling

Mutations and gametes go extinct during the course of a simulation.  For a gamete, this means that its count (`KTfwd::gamete_base::n`) is zero.  For a mutation, this means that the vector tracking mutation occurences has a zero at a specific position.

Internally, __fwdpp__ records where extinct mutations/gametes are and builds a FIFO queue of their indexes.  This queue is used so that these locations can be recycled with new objects.

Recycling allows us to use cache-friendly containers like `std::vector`.  It also allows us to re-use space allocated by gametes for their mutation keys.  Overall, it is a big win in terms of performance.

## Modeling the biology

This section gives a quick overview of what mutation, recombination, and genetic value functions have to look like in order to be compatible with __fwdpp__.  After reading this document, you should go look at the source code for the example programs to get an idea of how these ideas are put into practice.

This section uses the following shorthand notation for types:

* diploid_t is a diploid
* gamete_t is a gamete
* mutation_t is a mutation
* gcont_t is a container of gametes
* mcont_t is a container of mutations

In practice, the above will refer equally to both type names that are required to instantiate the library's template functions and to specific types used to instantiate the templates.  The following code block shows how `diploid_t` can be both a generic name and a specific type (`pair<size_t,size_t>`):

```cpp
template<typename diploid_t>
void do_something(const diploid_t &)
{
    static_assert(KTfwd::is_diploid<diploid_t>::value,
    "diploid_t must be a valid diploid type");
}

std::pair<std::size_t,std::size_t> a_diploid;

do_something(a_diploid);
```

The next sections discuss functions _signatures_.  By this, I mean its return type and the argument types that it expects. Further, the signature refers to the argument types _after_ lambda capture and/or binding parameters via `std::bind`.

The signature of our `do_something` function above can be written as `void(const diploid_t &)`.  I will refer to signatures using their representation as a `std::function`.  With this notation, we can rewrite the signature of `do_something` as

```cpp
std::function<void(const diploid_t &)>
```

In other words, anything convertible to the above type has the same signature.  In C++, this means regular functions, function objects defining `operator()`, lamba expressions, and C-style function pointers.

Consider the following function object:

```cpp
struct do_something_else
{
    template<typename diploid_t>
    void operator()(const diploid_t & dip, 
                    const double x) const
    {
    }
};
```

Further, consider the following object:

```cpp
auto bound_function = 
    std::bind(do_something_else(),
              std::placeholders::_1,2.0);
```

The signature of `bound_function` is  `std::function<void(const diploid_t &)>`. (In fact, that is one of many possible signatures, depending on what that placeholder eventually resolves to at compile time, but that detail is beyond the scope for now.)

Likewise, the following lambda has `std::function<void(const diploid_t &)>` for a signature:

```cpp
double x = 2.0;

auto a_lambda = 
    [x](const diploid_t &d) -> void
    {
    };
```

Much of __fwdpp__'s flexibility comes from the fact that many different callable objects can reduce to the same signature after binding/lambda capture.

### Mutation

Two types of mutation function are possible.  First,

```cpp
std::function<std::size_t(const recycling_bin_t &, mcont_t &)>
```

Or,

```cpp
std::function<std::size_t(const recycling_bin_t &, const gamete_t &,mcont_t &)>
```

The function must generate a new mutation and return its key.  The key must be the index of the new mutation in the mutations container.  Further, fwdpp's recycling rules must be obeyed.

In order to assist writing mutation functions, the library provids KTfwd::fwdpp_internal::recycle_mutation_helper.  This function is a variadic template function with the following prototype:

```cpp
template <typename queue_t, typename mcont_t, class... Args>
typename queue_t::value_type
recycle_mutation_helper(queue_t &mutation_recycling_bin,
                        mcont_t &mutations, Args &&... args)
```

The parameter pack Args represent the constructor calls for the mutation type, which are perfectly-forwarded into a new mutation object. If an existing location in the mutation container cannot be recycled, a new object is emplaced at the end using the parameter pack, minimizing the number of temporary objects created.  If a position can be recycled, a new object is most likely copy-constructed into place (but some compilers may implement a move construction at their discretion--most don't, though).

For a real-world example, see the implementation of KTfwd::infsites.

#### Possible future changes

* Move recycle_mutation_helper out of the fwdpp_internal namespace, as the public API should not be there.
* Consider additional valid function signatures where a diploid is passed in, which may help implement things like sex-specific mutation rates, etc.
* Design an API that allows for reversible mutation.  I'm not sure fwdpp is quite there yet.

### Recombination

A recombination function has the following signature:

```cpp
std::function<std::vector<double>(const gamete_t &, const gamete_t &, const mcont_t &)>
```

This function is responsible for returning a vector of breakpoints.  Any random number generators, etc., needed to generate such a vector must be bound/captured separately.

An example of a recombination function is KTfwd::poisson_xover, which is implemented as a function object with a template call operator.

This function signature is the same for a multilocus diploid.  Vectors of such functions are used for multi-locus/region simulations, and they are applied in turn to each locus.

#### Possible future changes

In order to model things like sex-limited recombination, etc., will require a signature like this:

```cpp
std::function<std::vector<double>(const diploid_t &, const gcont_t &, const mcont_t &)>
```

Passing in the diploid type itself would allow checking of any data present in custom diploid types (male, female, etc.).  It is likely that such a change will happen in a future release.

### Calculating a diploid's genetic value

For simulations with fitness effects, we need to:

1. Take the data in a diploid's gametes and calculate some quantity.  Call it the genetic value, which is a `double`.
2. Convert that value into fitness, which is a non-negative `double`.

For standard population-genetic simulations, the mapping of genetic value to fitness is trivial.  Often, it is as simple as `max(0.0,genetic_value);`.

For simulations of quantitative traits, there are probably extra steps involved:

1. Is the final trait value going to be a combination of genetic value plus random noise?
2. What is the mapping of final trait value to fitness?  For example--Gaussian stabilizing selection or a linear selection gradient?

For quantitative genetic simulations, the genetic value is a `double` that can take on any finite value.

We define the following conventions:

1. Fitnesses must be `w >= 0` and a mutant-free diploid would have a fitness of `w = 1`. 
2. Genetic values may take on any finite value and a mutant-free diploid would have `g = 0`.

The challenge is to write functions to calculate `g` and map them onto the correct scale.

__fwdpp__ has a simple requirement for the signature of a function to calculate a genetic value:

```cpp
std::function<double(const diploid_t &, const gcont_t &, const mcont_t &)>
```

In other words, it take a diploid, a container of gametes, and a container of mutations as arguments.  Using that information, the genetic value is calculated and returned as a double.

For the case of a single locus/region simulation, the library provides efficient implementations of two standard models:

* `KTfwd::additive_diploid` applies an additive model with dominance.
* `KTfwd::multiplicative_diploid` applies a multiplicative model with dominance.

By default, both of these function objects map genetic values to fitness when constructing the return value.  This is done via calls to `KTfwd::aw` or `KTfwd::mw`, respectively.  This behavior may be changed by passing in a differnt policy to the object's constructor--`KTfwd::atrait()` and `KTfwd::mtrait()`, respectively, change the behavior of each type to return a trait value centered on zero.

#### General comments

Mapping a diploid genotype to fitness is one of the most important aspects of implementing a simulation.  It is also one of the most difficult, both in terms of the programming itself and because it really forces you to think about your models.  For "standard pop-gen" scenarios, the default behaviors are quite simple.  For quantitative trait simulations, things are more necessarily more complex, and the library allows you to do any of the following:

1. Map genetic value straight to fitness.
2. Separate genetic value from trait value from fitness calculations.

The latter is very important, as it allows the flexibility of recording the `g`, `e`, and `w` values that may be part of the calculation as data, for example in a custom diploid type.

For more details, see the tutorials.
