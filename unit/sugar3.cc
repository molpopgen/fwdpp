/*!
  \file sugar3.cc
  \ingroup unit
  \brief Testing KTfwd::multiloc and KTfwd::multiloc_serialized
*/
#define BOOST_TEST_MODULE sugarTest3
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <functional>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#define FWDPP_SUGAR_USE_BOOST
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;

//Fitness function
struct no_selection_multi
{
  typedef double result_type;
  template< typename diploid_type >
  //inline double operator()( const std::vector< std::pair<iterator_type,iterator_type> > & diploid ) const
  inline double operator()( diploid_type & diploid ) const
  {
    return 1.;
  }
};

BOOST_AUTO_TEST_CASE( metapop_sugar_test1 )
{
  using poptype = KTfwd::multiloc<KTfwd::popgenmut>;
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  poptype pop(1000,4);
  unsigned generation = 0;
  
  //Mutation model for 4 loci
  std::vector< std::function<KTfwd::popgenmut(typename poptype::mlist_t *)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,generation,
	      0.005,0.,[](gsl_rng * r){return gsl_rng_uniform(r);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,1.,2.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,2.,3.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),rng,std::placeholders::_1,&pop.mut_lookup,generation,
      		0.005,0.,[](gsl_rng * r){return gsl_ran_flat(r,3.,4.);},[](gsl_rng * r){return 0.;},[](gsl_rng * r){return 0.;})
  };
  
  //Within-locus recombination models for 4 loci
  std::vector< std::function<unsigned(typename poptype::glist_t::iterator &,
				      typename poptype::glist_t::iterator &)> > recmodels {
    //Locus 0: positions = uniform [0,1)
    std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[0],0.005,rng,[&rng](){ return gsl_rng_uniform(rng); }),
      //Locus 1: positions = uniform [1,2)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[1],0.005,rng,[&rng](){ return gsl_ran_flat(rng,1.,2.); }),
      //Locus 2: positions = beta(1,10) from 2 to 3
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[2],0.005,rng,[&rng](){ return 2. + gsl_ran_beta(rng,1.,10.); }),
      //Locus 2: positions = uniform [3,4)
      std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&pop.gametes[3],0.005,rng,[&rng](){ return gsl_ran_flat(rng,3.,4.); })
      };

  //Equal mutation and rec. rates per locus
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng,
					   &pop.gametes,
					   &pop.diploids,
					   &pop.mutations,
					   1000,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(KTfwd::insert_at_end<poptype::mutation_t,poptype::mlist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(KTfwd::insert_at_end<poptype::gamete_t,poptype::glist_t>,std::placeholders::_1,std::placeholders::_2),
					   std::bind(no_selection_multi(),std::placeholders::_1),
					   std::bind(KTfwd::mutation_remover(),std::placeholders::_1,0,2000),
					   0.);
      for( unsigned i = 0 ; i < 4 ; ++i )
	{
	  assert( check_sum(pop.gametes[i],2000) );
	}
      KTfwd::remove_fixed_lost(&pop.mutations,&pop.fixations,&pop.fixation_times,&pop.mut_lookup,generation,2000);
    }
}
