# Tutorial 5: the "sugar" layer on top of fwdpp

Through __fwdpp__ 0.2.9, a lot of the basic setup of a simulation could get tedious due to having to define a large number of container types.  In version 0.3.0 of the library, I introduced the subdirectory fwdpp/sugar, which is a layer of "syntatic sugar" designed to make fwdpp both easier to use and easier to wrap for use in other programming environments.

These functions are only relevant to _indivdual-based_ simulations.  The gamete-based portion of __fwdpp__ is unlikely to be supported (which is fine, as it is slower...).

## Purpose

The intention of the "sugar layer" is to:

* Provide a path to more rapid development of simulations, which is accomplished by providing a set of C++11 template aliases for structures defining population objects.
* Provide a set of standard mutation types and a standard mutation model policy to support such objects.  These types include fwdpp::mutation (which has been in the library "forever"), fwdpp::popgenmut, which is implemented overm, which is used in the examples (@ref md_md_examples), and fwdpp::infsites_popgenmut, which implements an infinitely-many sites mutation scheme generate the latter mutation type. 
* Provide data types wrapping the serialization methods in __fwdpp__ (see @ref md_md_serialization).  The relevant sugar components are fwdpp::serialize and fwdpp::deserialize.
* Make sure that all types defined in the sugar layer are properly-defined so that they may be used in conjunction with tools like [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) to enable running/processing simulations within R, [boost](http://www.boost.org)'s "boost.python" library to enable processing/running simulations in python, and [pybind11](https://github.com/pybind/pybind11), a C++11/c++14 reimagining of boost.python.

## Overview

The sugar layer is comprised of several different headers, all located in _fwdpp/sugar_.  These header files are intended to be used directly for implementing simulations.

Further subdirectories contain the details of the namespace fwdpp::sugar, which most library users are unlikely to need to deal with directly.

Importantly, _the use of these headers is optional_, and may not work for every type of simulation.  Much of this is based upon what I've found to work for my own research using __fwdpp__.

## Components 
### Smart wrapper to GSL random number generators

__fwdpp__ using the [GSL](http://gnu.org/software/gsl) for random number generation.  The sugar layer provides a smart pointer wrapper around a gsl_rng *.  Currently, the only supported rng types are gsl_rng_mt19937 and gsl_rng_taus2.

The wrapper type is fwdpp::GSLrng_t:

~~~{.cpp}
#include <fwdpp/sugar/GSLrng_t.hpp>
~~~

To specify the use of the Mersenne twister:

~~~{.cpp}
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;
GSLrng rng(101); //101 is the seed value 
~~~

The taus generator may be declared like this:

~~~{.cpp}
using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2>;
GSLrng rng(101); //101 is the seed value 
~~~

This wrapper has the following properties:
* It is constructed with an unsigned integer, which will be used as a seed
* It is implemented in terms of a std::unique_ptr
* It is not copyable
* The underlying pointer is retrived via a call to get()
* It frees the gsl_rng * upon destruction using fwdpp::sugar::gsl_rng_deleter

The template parameters are of type fwdpp::sugar::GSL_RNG_TYPE_TAG;

### Built-in "standard" mutation types and models.

__fwdpp__ provides fwdpp::mutation, which is a "standard" type of variant for a population-genetic simulation.  A mutation is association with a selection coefficient, \f$s\f$, and a dominance coefficient, \f$h\f$.  The typical "popgen" type of simulation assigns traits values \f$1\f$, \f$1+hs\f$, and \f$1+xs\f$ (where \f$x = 1\f$ or 2, depending on the particulars) to genotypes \f$AA\f$, \f$Aa\f$, and \f$aa\f$, respectively, where \f$a\f$ represents the mutation.

The sugar layer provides the additional type fwdpp::popgenmut, which also tracks \f$g\f$, the generation in the simulation where the mutation first arose.

The relevant sugar headers are:

~~~{.cpp}
//for fwdpp::popgenmut
#include <fwdpp/sugar/popgenmut.hpp>
~~~

### Simplifying the declaration of a population

The sugar layer greatly simplifies the setup of the various containers needed for a simulation.  The typical simulation requires that the programmer define mutation types, gamete types, mutation container types, gamete container types, and diploid vector types.  Further, the programmer may choose to use boost's rapid memory allocation pools and use conditional compilation to support systems without boost.

The sugar layer allows you to set up a simulation only having to decide on your mutation type.  The types of simulations that are supported are given by the enum type fwdpp::sugar::FWDPP_SUGAR_POPTYPE and the dispatch tag type fwdpp::sugar::FWDPP_SUGAR_POPTAG.

Let's assume that you want to do a simulation implemented in terms of fwdpp::popgenmut.  You may set up various types of simulations like this:

~~~{.cpp}
//Single-locus population
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>

using poptype = fwdpp::slocuspop<fwdpp::popgenmut>;
~~~


~~~{.cpp}
//Multi-locus population
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/mlocuspop.hpp>

using poptype = fwdpp::mlocuspop<fwdpp::popgenmut>;
~~~

The code chunks above will result in the following:

* All containers are from namespace std
* Populations are copy- and move-constructible
* Populations are serializable via fwdpp::io::serialize_population

#### Details

These "poptypes" contain all of the types needed to make calls to fwdpp::sample_diploid, etc.  

See the following files for concrete working examples:

* diploid_ind.cc
* diploid_fixed_sh_ind.cc

#### Further customization

It you aren't happy with how I've set up the "poptypes", you may provide your own typedefs in terms of the following classes:

* fwdpp::sugar::slocuspop
* fwdpp::sugar::mlocuspop 

You are able to publicly inherit them or encapsulate them in the usual ways, if more customization is needed.

__fwdpp__ 0.3.1 added support for custom diploid types (see @ref md_md_customdip).  The template aliases for population types were updated accordingly.

