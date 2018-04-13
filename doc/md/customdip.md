# Advanced topics: custom diploid genotypes

## Intro

This section discusses how to implement simulations where a diploid's genotype is represented by a user-defined type.

This is a document covering advanced programming issues using __fwdpp__.  Please see the introductory tutorials if you are new to programming with the library.

## Rationale

A minimal diploid contains two integers representing the locations of the diploid's gametes.  These types should be the
same as the size_type corresponding to the gamete container.  Typically, std::size_t suffices, and thus a minimal
diploid may be defined as:

~~~{.cpp}
using diploid_t = std::pair<std::size_t,std::size_t>;
~~~

It is often desirable to have diploid types with additional data, and __fwdpp__ easily supports this in the form of "custom" diploid types. At a minimum, a custom diploid must have the same public API as 
[std::pair](http://en.cppreference.com/w/cpp/utility/pair). In other words, it must look like this at a minimum:

~~~{.cpp}
struct my_custom_diploid
{
    using first_type = std::size_t;
    using second_type = std::size_t;
    first_type first;
    second_type second;
};
~~~

To check that a custom diploid type is acceptable, you can use fwdpp's type traits library:

~~~{.cpp}
#include <fwdpp/traits.hpp>

static_assert(fwdpp::traits::is_diploid<my_custom_diploid>::value,"Oops!");
~~~

An additional type trait function exists to detect custom diploids (fwdpp::traits::is_custom_diploid).  It checks that the minimum API requirements are met, but that the type is not identical to a std::pair.

For a fully worked-out example, see @ref custom_diploid.cc
