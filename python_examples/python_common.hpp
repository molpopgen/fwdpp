#ifndef __FWDPP_EXTENSIONS_PYTHON_COMMON_HPP__
#define __FWDPP_EXTENSIONS_PYTHON_COMMON_HPP__

//include main fwdpp library
#include <fwdpp/diploid.hh>
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

using poptype = KTfwd::singlepop<KTfwd::mutation>;

//Evolve the population for some amount of time with mutation and recombination
poptype evolve( GSLrng & rng,
		const unsigned & N,
		const unsigned & generations,
		const double & mu,
		const double & recrate)
{
  poptype pop(N);
  pop.mutations.reserve(std::ceil(std::log(2*N)*(4.*double(N)*(mu))+0.667*(4.*double(N)*(mu))));
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get()); //uniform crossover map
  for( unsigned generation = 0 ; generation < generations ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng.get(),
					  pop.gametes,
					  pop.diploids,
					  pop.mutations,
					  pop.mcounts,
					  N,
					  mu,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),
						    mu,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
					  std::bind(KTfwd::poisson_xover(),rng.get(),recrate,0.,1.,
						    std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
						    std::placeholders::_3,2.),
					  pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*pop.N);
    }
  return pop;
}

#endif
