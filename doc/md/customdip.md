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

It is often desirable to have diploid types with additional data, and __fwdpp__ easily supports this in the form of
"custom" diploid types.  Such custom types must have the same public API as
[std::pair](http://en.cppreference.com/w/cpp/utility/pair). Further, the type must inherit from
KTfwd::tags::custom_diploid_t.  Once these requirements are fulfilled, the type will be
compatible with the library. See type_traitsTest.cc for minimal examples.
