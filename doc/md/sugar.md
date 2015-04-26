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

The sugar layer also provides KTfwd::infsites, which generates mutations according to an infinitely-many sites scheme.  The mutation types that are supported are KTfwd::mutation and KTfwd::popgenmut.  I intend KTfwd::infsites to facilitate rapid prototyping of __fwdpp__-based simulations, and many users may find these mutation-related types sufficient for their own research needs.

The relevant sugar headers are:

~~~{.cpp}
//for KTfwd::popgenmut
#include <fwdpp/sugar/popgenmut.hpp>
//for KTfwd::infsites
#include <fwdpp/sugar/infsites.hpp>
~~~

## Simplifying the declaration of a population

The sugar layer greatly simplifies the setup of the various containers needed for a simulation.  The typical simulation requires that the programmer define mutation types, gamete types, mutation list types, gamete list types, and diploid vector types.  Further, the programmer may choose to use boost's rapid memory allocation pools and use conditional compilation to support systems without boost.

The sugar layer allows you to set up a simulation only having to decide on your mutation type.  The types of simulations that are supported are given by the enum type KTfwd::sugar::FWDPP_SUGAR_POPTYPE and the dispatch tag type KTfwd::sugar::FWDPP_SUGAR_POPTAG.

Let's assume that you want to do a simulation implemented in terms of KTfwd::popgenmut.  You may set up various types of simulations like this:

~~~{.cpp}
//Single population, single locus
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>

using poptype = KTfwd::singlepop<KTfwd::popgenmut>;
~~~

~~~{.cpp}
//Metapopulation, single locus
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/metapop.hpp>

using poptype = KTfwd::metapop<KTfwd::popgenmut>;
~~~

~~~{.cpp}
//Single population, multi locus
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/multiloc.hpp>

using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
~~~

The code chunks above will result in the following:

* All containers are from namespace std
* Populations are move-constructible
* Populations are __not__ copyable.  The copy constructor and assignment operator are both deleted.

In order to use boost containers and boost::pool_allocator instead of the std containers, define FWDPP_SUGAR_USE_BOOST for your preprocessor.  For example:

~~~{.cpp}
//Single population, single locus
#include <fwdpp/sugar/popgenmut.hpp>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/singlepop.hpp>

using poptype = KTfwd::singlepop<KTfwd::popgenmut>;
~~~

(Note that simulations will typically be a lot faster if you use boost's memory pools!)

### Copyable populations

There are important use cases that require copyable populations.  For example, exposing classes to python using boost.python requires that the underlying C++ types are copy-constructible, and thus the "poptypes" described above will not compile successfully because their copy constructors are deleted for safety.

However, due to the pointer/iterator-based design of __fwdpp__, such copies require deep copies that reconstruct data structures.  Such copies are complex, and can be handled manually via KTfwd::write_binary_pop and KTfwd::write_binary_metapop to _write_ data to an internal buffer (such as std::ostringstream), and KTfwd::read_binary_pop and KTfwd::read_binary_metapop to _read_ populations from a buffer.  These same functions can be used to implement copy constructors and assignment operators for classes representing populations.  However, such objects require policies to write and read the simulation's mutation type (see \ref md_md_serialization for how to write such policies for custom mutation types).

Once such serialization policies are implemented, you may use the following types that _are_ copy-constructible and assignable via operator=:

* KTfwd::singlepop_serialized
* KTfwd::metapop_serialized
* KTfwd::multiloc_serialized

Now, copy construction will work:

~~~{.cpp}
#include <fwdpp/sugar/serialization.hpp>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/singlepop.hpp>

using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
using poptype = KTfwd::singlepop_serialized<KTfwd::popgenmut,mwriter,mreader>;
int main( int argc, char ** argv )
{
	//N = 100
	poptype pop(100);
	poptype pop2(pop);
}
~~~

### Details

These "poptypes" contain all of the types needed to make calls to KTfwd::sample_diploid, etc.:

~~~{.cpp}
#include <fwdpp/diploid.hh>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mutation_with_age = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_with_age>;
using poptype = KTfwd::singlepop<mutation_with_age>;

int main(int argc, char ** argv)
{
  poptype pop(1000);

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng);
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng,
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,generation,
						    0.005,0.,[](gsl_rng * r){return gsl_rng_uniform(r);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						    &pop.gametes,
						    0., //no rec
						    rng,
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop.N));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
  }
~~~

## Simplifying serializing of simulated data

The population types discussed above may be serialized using KTfwd::serialize and deserialized using KTfwd::deserialize.  The sugar layer also provides KTfwd::mutation_writer and KTfwd::mutation reader to support the serialization of KTfwd::mutation and KTfwd::popgenmut.

For example:

~~~{.cpp}
#include <fwdpp/sugar/serialization.hpp>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/singlepop.hpp>

using poptype = KTfwd::singlepop<KTfwd::popgenmut>;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
int main( int argc, char ** argv )
{
	//N = 100
	poptype pop(100);
	KTfwd::serialize s;
	s(pop,mwriter);
	poptype pop2(0);
	//Copy pop into pop2 via the buffer in s
	KTfwd::deserialize()(pop2,s,mreader());
}
~~~

The above will also work with the _serialized versions of the various "poptypes".

### Details

KTfwd::serialize contains a std::stringstream called "buffer", which will hold the serialized data.  After deserializing, buffer.seekg(0) is called so that you may read from the beginning of the buffer.

You may also use buffer to write to a file:

~~~{.cpp}
#include <fwdpp/sugar/serialization.hpp>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/singlepop.hpp>
#include <zlib.h>
using poptype = KTfwd::singlepop<KTfwd::popgenmut>;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;
int main( int argc, char ** argv )
{
	//N = 100
	poptype pop(100);
	KTfwd::serialize s;
	s(pop,mwriter);
	gzFile gzf = gzopen("pop.gz","wb");
	gzwrite( gzf, s.buffer.str(), s.buffer.str().size() );
	gzclose(gzf);
}
~~~

You may use KTfwd::serialize and KTfwd::deserialize for any custom mutation types, once the proper serializing objects are written (see \ref md_md_serialization).

In future releases, support for direct serialization to/from gzipped files will be added.  In the mean time, you can do this manually using __fwdpp__ functions like KTfwd::write_binary_pop, etc.
