# Preprocessor symbols

This document describes the various preprocessor symbols affecting how simulations compiled against __fwdpp__ behave.

## Backwards compatibility/reproducibility

With the exception of bug fixes, the intention is that simulation results are reproducible when compiled against different versions of the library.  Sometimes, though, algorithm improvements may be introduced that change the order and/or the number of calls to a random number generator.  Such changes necessarily change the final output of a simulation.  Where possible, the preprocessor can be used to select between the new/improved code and the older code for situations where you want to reproduce results exactly (assuming the same GSL random number generator, seed, etc.).

The format for these symbols is COMPAT_X, meaning that defining the symbols results in code identical to that used in __fwdpp__ versions _less than or equal to X_.

The symbols curently in use are:

* COMPAT_0_3_0: version 0.3.1 of the library made changes to the implementation of selfing.  If you need to reproduce results based on __fwdpp__ version prior to 0.3.1, compile your program with -DCOMPAT_0_3_0

### Example use

To compile the example programs using the selfing implementation from version 0.3.0 and earlier, you may configure the library like this:

~~~{.sh}
./configure CXXFLAGS="-O2 -g -DCOMPAT_0_3_0"
~~~

You need to add the "-g -O2" as well, because you are asking the configure script to override the default compiler flags with your own custom flags, but you definitely want to compile with optimization!

## Internal use of boost or std containers

Testing during development identified some cases where there is a performance advantage to using boost::vector instead of std::vector.  The decision of what to use is controlled by code like this:

~~~{.cpp}
#if defined(HAVE_BOOST_VECTOR) && !defined(USE_STANDARD_CONTAINERS)
//Use boost::container::vector
#else
//use std::vector
#endif
~~~

The most straightforward/portable way to deal with these two symbols is via a system like [GNU autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html).  Specifically, you want the ./configure step to output a config.h header that your programs will include.  The easiest way to do this is to use the "devtools" included with the library source code when starting your own projects based on __fwdpp__  (@ref md_md_devtools).
