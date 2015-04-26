# Tutorial 5: the "sugar" layer on top of fwdpp

Through __fwdpp__ 0.2.9, a lot of the basic setup of a simulation could get tedious due to having to define a large number of container types.  In version 0.3.0 of the library, I introduced the subdirectory fwdpp/sugar, which is a layer of "syntatic sugar" designed to make fwdpp both easier to use and easier to wrap for use in other programming environments.

These functions are only relevant to _indivdual-based_ simulations.  The gamete-based portion of __fwdpp__ is unlikely to be supported (which is fine, as it is slower...).

## Overview

The sugar layer is comprised of several different headers, all located in _fwdpp/sugar_.  These header files are intended to be used directly for implementing simulations.

Further subdirectories contain the details of the namespace KTfwd::sugar, which most library users are unlikely to need to deal with directly.

Importantly, _the use of these headers is optional_, and may not work for every type of simulation.  Much of this is based upon what I've found to work for my own research using __fwdpp__.

A major goal of this part of __fwdpp__ is to facilitate using the library in other programming environments, especially R and python.  I have successfully used this code to execute simulations in R using [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) and in python via [boost](http://www.boost.org)'s "boost.python" library.

You may find an example of using __fwdpp__ + Rcpp in [foRward](http://github.com/molpopgen/foRward) and an example of using boost.python in the extensions subdirectory of the __fwdpp__ source code repository.


## Smart wrapper to GSL random number generators

__fwdpp__ using the [GSL](http://gnu.org/software/gsl) for random number generation.  The sugar layer provides a smart pointer wrapper around a gsl_rng *.  Currently, the only supported rng types are gsl_rng_mt19937 and gsl_rng_taus2.

The wrapper type is KTfwd::GSLrng_t:

~~~{.cpp}
#include <fwdpp/sugar/GSLrng_t.hpp>
~~~

To specify the use of the Mersenne twister:

~~~{.cpp}
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;
GSLrng rng(101); //101 is the seed value 
~~~

The taus generator may be declared like this:

~~~{.cpp}
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2>;
GSLrng rng(101); //101 is the seed value 
~~~

This wrapper has the following properties:
* It is constructed with an unsigned integer, which will be used as a seed
* It is implemented in terms of a std::unique_ptr
* It is copy-constructable via a call to gsl_rng_clone to create a new unique_ptr
* It is implicitly convertible to a gsl_rng *, and thus may be passed to any function taking such a type as a parameter
* It frees the gsl_rng * upon destruction using KTfwd::sugar::gsl_rng_deleter

The template parameters are of type KTfwd::sugar::GSL_RNG_TYPE_TAG;

## Built-in "standard" mutation types and models.

__fwdpp__ provides KTfwd::mutation, which is a "standard" type of variant for a population-genetic simulation.  A mutation is association with a selection coefficient, $s$, and a dominance coefficient, $h$.  The typical "popgen" type of simulation assigns traits values $1$, $1+hs$, and $1+xs$ (where $x = 1$ or 2, depending on the particulars) to genotypes $AA$, $Aa$, and $aa$, respectively, where $a$ represents the mutation.

The sugar layer provides the additional type KTfwd::popgenmut, which also tracks $g$, the generation in the simulation where the mutation first arose.

## Simplifying the declaration of a population

## Simplifying serializing of simulated data


