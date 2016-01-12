/*! \file sugar2.cc
  \ingroup unit 
  \brief Testing KTfwd::metapop and KTfwd::metapop_serialized
*/
#define BOOST_TEST_MODULE sugarTest2
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mutation_t = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_t>;

size_t migpop(const size_t & source_pop, gsl_rng * r, const double & mig_prob)
{
  if( gsl_rng_uniform(r) < mig_prob )
    {
      return ! source_pop;
    }
  return source_pop;
}

using poptype = KTfwd::metapop<mutation_t>;

void simulate(poptype & pop)
{
  //Evolve for 10 generations
  std::vector<std::function<double (const poptype::gamete_t &,
				    const poptype::gamete_t &,
				    const poptype::mcont_t &)> > fitness_funcs(2,
									       std::bind(KTfwd::multiplicative_diploid(),std::placeholders::_1,std::placeholders::_2,
											 std::placeholders::_3,2.));
  KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
  for( unsigned generation= 0 ; generation < 10 ; ++generation )
    {
      std::vector<double> wbar = KTfwd::sample_diploid(rng.get(),
						       pop.gametes,
						       pop.diploids,
						       pop.mutations,
						       pop.mcounts,
						       &pop.Ns[0],
						       0.005,
						       std::bind(KTfwd::infsites(),std::placeholders::_1,std::placeholders::_2,rng.get(),std::ref(pop.mut_lookup),generation,
								 0.005,0.,[&rng](){return gsl_rng_uniform(rng.get());},[](){return 0.;},[](){return 0.;}),
						       std::bind(KTfwd::poisson_xover(),rng.get(),0.005,0.,1.,
								 std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
						       fitness_funcs,
						       std::bind(migpop,std::placeholders::_1,rng.get(),0.001),
						       pop.neutral,pop.selected);
      KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,generation,4000);
    }
}

BOOST_AUTO_TEST_CASE( metapop_sugar_test1 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2(pop);
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( metapop_sugar_test2 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2{0,0};
  KTfwd::serialize s;
  s(pop,mwriter());
  KTfwd::deserialize()(pop2,s,mreader());
  BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE( metapop_sugar_test3 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2(std::move(pop));
  BOOST_CHECK_EQUAL(pop==pop2,false);
}

BOOST_AUTO_TEST_CASE( metapop_sugar_test4 )
{
  poptype pop({1000,1000});
  simulate(pop);
  poptype pop2=std::move(pop);
  BOOST_CHECK_EQUAL(pop==pop2,false);
}
