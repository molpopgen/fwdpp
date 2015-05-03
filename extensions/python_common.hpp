#ifndef __FWDPP_EXTENSIONS_PYTHON_COMMON_HPP__
#define __FWDPP_EXTENSIONS_PYTHON_COMMON_HPP__

//include main fwdpp library
#include <fwdpp/diploid.hh>
//Tell the fwdpp "sugar" layer to use boost containers (no harm here, as it is likely installed b/c we need boost python..)
#define FWDPP_SUGAR_USE_BOOST
//Include the necessary "sugar" components
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/infsites.hpp>
/*
  We will use a gsl_rng_mt19937 as our RNG.
  This type is implicitly convertible to gsl_rng *,
  and auto-handles the gsl_rng_free steps, etc.
 */
using GSLrng = KTfwd::GSLrng_t<KTfwd::GSL_RNG_MT19937>;

/*
  boost.python requires that any C++ classes exposed to python
  are copy-constructible.

  As far as fwdpp is concerned, this requirement means that 
  we must use the serializable version of the fwdpp sugar libary.

  This example is implemented in terms of KTfwd::mutation, and
  the fwdpp sugar layer supports serialization of this type
 */
using poptype = KTfwd::singlepop_serialized<KTfwd::mutation,KTfwd::mutation_writer,KTfwd::mutation_reader<KTfwd::mutation>>;

//Evolve the population for some amount of time with mutation and recombination
poptype evolve( GSLrng & rng,
		const unsigned & N,
		const unsigned & generations,
		const double & mu,
		const double & recrate)
{
  poptype pop(N);
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng); //uniform crossover map
  for( unsigned generation = 0 ; generation < generations ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng,
					  &pop.gametes,
					  &pop.diploids,
					  &pop.mutations,
					  N,
					  mu,
					  std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,
						    mu,0.,[&]{return gsl_rng_uniform(rng);},[]{return 0.;},[]{return 0.;}),
					  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,
						    &pop.gametes,
						    recrate, 
						    rng,
						    recmap),
					  std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,2.),
					  std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2*pop.N));
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2*pop.N);
    }
  return pop;
}

#endif
