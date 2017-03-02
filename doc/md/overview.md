# Overview of fwdpp features

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
