/*!
  \file sugar_multilocus.cc
  \ingroup unit
  \brief Testing KTfwd::multiloc
*/

#define BOOST_TEST_MODULE sugar_multilocus
#define BOOST_TEST_DYN_LINK

#include <unistd.h>
#include <config.h>
#include <zlib.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>

using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<KTfwd::popgenmut>;

using poptype = KTfwd::multiloc<KTfwd::popgenmut>;

//Fitness function
struct multilocus_additive
{
  using result_type = double;
  inline double operator()(const poptype::dipvector_t::value_type & diploid,
			   const poptype::gcont_t & gametes,
			   const poptype::mcont_t & mutations) const
  {
    using dip_t = poptype::dipvector_t::value_type::value_type;
    return std::max(0.,1.+ std::accumulate( diploid.begin(),diploid.end(),0.,[&gametes,&mutations](const double d,const dip_t & dip)
				{
				  return d + KTfwd::additive_diploid()(gametes[dip.first],gametes[dip.second],mutations) - 1.;
				} ));
  }
};

void simulate( poptype & pop )
{

  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  unsigned generation = 0;

  std::vector<std::function<std::size_t(std::queue<std::size_t> &,poptype::mcont_t &)> > mutmodels {
    //Locus 0: positions Uniform [0,1)
    std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),&generation,
	      0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}) ,
      //Locus 1: positions Uniform [1,2)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),1.,2.);},[](){return 0.;},[](){return 0.;}),
      //Locus 2: positions Uniform [2,3)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),2.,3.);},[](){return 0.;},[](){return 0.;}),
      //Locus 3: positions Uniform [3,4)
      std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),&generation,
      		0.005,0.,[&rng](){return gsl_ran_flat(rng.get(),3.,4.);},[](){return 0.;},[](){return 0.;})
      };

  //Within-locus recombination models for 4 loci
  std::vector<std::function<std::vector<double>(const poptype::gamete_t &,
						const poptype::gamete_t &,
						const poptype::mcont_t &)> > recmodels
  {
    std::bind(KTfwd::poisson_xover(),rng.get(),0.005,0.,1.,
	      std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
      std::bind(KTfwd::poisson_xover(),rng.get(),0.005,1.,2.,
		std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
      std::bind(KTfwd::poisson_xover(),rng.get(),0.005,2.,3.,
		std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
      std::bind(KTfwd::poisson_xover(),rng.get(),0.005,3.,4.,
		std::placeholders::_1,std::placeholders::_2,std::placeholders::_3)
      };

  //Equal mutation and rec. rates per locus
  std::vector<double> mu(4,0.005),rbw(3,0.005);
  BOOST_REQUIRE_EQUAL(mu.size(),4);
  for( ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid( rng.get(),
					   pop.gametes,
					   pop.diploids,
					   pop.mutations,
					   pop.mcounts,
					   1000,
					   &mu[0],
					   mutmodels,
					   recmodels,
					   &rbw[0],
					   [](gsl_rng * __r, const double __d){ return gsl_ran_binomial(__r,__d,1); },
					   std::bind(multilocus_additive(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					   pop.neutral,pop.selected);
      assert( check_sum(pop.gametes,8000) );
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2000);
    }
}

BOOST_AUTO_TEST_CASE( multiloc_sugar_test1 )
{
  poptype pop(1000,4);
  simulate(pop);
  poptype pop2(pop);
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( multiloc_sugar_test2 )
{
  poptype pop(1000,4);
  simulate(pop);
  poptype pop2(0,0);
  KTfwd::serialize s;
  s(pop,mwriter());
  KTfwd::deserialize()(pop2,s,mreader());
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( multiloc_sugar_test2_gz )
{
  poptype pop(1000,4);
  simulate(pop);
  poptype pop2(0,0);
  gzFile gzf = gzopen("sugar_multilocus_out.gz","wb");
  KTfwd::gzserialize()(gzf,pop,mwriter());
  gzclose(gzf);
  gzf=gzopen("sugar_multilocus_out.gz","rb");
  KTfwd::gzdeserialize()(gzf,pop2,std::bind(mreader(),std::placeholders::_1));
  gzclose(gzf);
  unlink("sugar_multilocus_out.gz");
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( multiloc_sugar_test3 )
{
  poptype pop(1000,4);
  simulate(pop);
  poptype pop2(std::move(pop));
  BOOST_CHECK_EQUAL(pop==pop2,false);
}

BOOST_AUTO_TEST_CASE( multiloc_sugar_test4 )
{
  poptype pop(1000,4);
  simulate(pop);
  poptype pop2=std::move(pop);
  BOOST_CHECK_EQUAL(pop==pop2,false);
}
