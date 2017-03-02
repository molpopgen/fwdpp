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

### Diploids

When modeling a single contiguous genomic segment, a diploid is simply a pair of gametes.  The way to represent that using __fwdpp__ is actually with a type from the C++ standard library:

```{.cpp}
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

Why does __fwdpp__ not simply have custom diploid types inherit from the standard pair type?  Well, the standard library pair does not have a virtual destructor, meaning that it cannot be used as a base type.  Future versions of the library may add a custom diploid base type.

#### Multi-locus/multi-region diploids

### Gametes

### Mutations

## Representing a population

## Modeling the biology

### Mutation

### Recombination

### Fitness
