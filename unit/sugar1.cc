/*! 
  \file sugar1.cc 
  \ingroup unit 
  \brief Testing KTfwd::singlepop 
*/
#define BOOST_TEST_MODULE sugarTest1
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>

using mutation_t = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_t>;
using singlepop_t = KTfwd::singlepop<mutation_t>;

void simulate( singlepop_t & pop )
{
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

  //Evolve for 10 generations
  std::function<double(void)> recmap = std::bind(gsl_rng_uniform,rng.get());
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      double wbar = KTfwd::sample_diploid(rng.get(),
					  pop.gametes,
					  pop.diploids,
					  pop.mutations,
					  pop.mcounts,
					  1000,
					  0.005,
					  std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),generation,
						    0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
					  std::bind(KTfwd::poisson_xover(),rng.get(),0.005,0.,1.,
						    std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
					  std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
						    std::placeholders::_3,2.),
					  pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,2*pop.N);
    }
}
  

BOOST_AUTO_TEST_CASE( singlepop_sugar_test1 )
{
  singlepop_t pop(1000);
  simulate(pop);

  singlepop_t pop2(pop);

  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_test2 )
{
  singlepop_t pop(1000);
  simulate(pop);

  KTfwd::serialize s;
  s(pop,mwriter());
  singlepop_t pop2(0);
  KTfwd::deserialize()(pop2,s,mreader());
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_test3 )
{
  singlepop_t pop(1000);
  simulate(pop);

  singlepop_t pop2(std::move(pop));
  //Should be false b/c move will leave
  //pop's containers in a wacky state
  BOOST_CHECK_EQUAL(pop==pop2,false);
}

BOOST_AUTO_TEST_CASE( singlepop_sugar_test4 )
{
  singlepop_t pop(1000);
  simulate(pop);

  singlepop_t pop2=std::move(pop);
  //Should be false b/c move will leave
  //pop's containers in a wacky state
  BOOST_CHECK_EQUAL(pop==pop2,false);
}
