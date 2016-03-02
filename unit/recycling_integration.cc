/*! 
  \file recycling_integration.cc 
  \ingroup unit 
  \brief Testing that recycling and cleanup functions work
*/
#define BOOST_TEST_MODULE sugar_sampling
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

using mutation_t = KTfwd::popgenmut;
using singlepop_t = KTfwd::singlepop<mutation_t>;

/*
  A note on test names:
  All recycling is applied to extinct mutations.
  What varies from test to test is the treatment of fixations.
*/

BOOST_AUTO_TEST_CASE( only_recycle_extinct_mutations )
{
  unsigned N=1000;
  singlepop_t pop(N);
  //This is a neutral mutation
  pop.mutations.emplace_back(0.1,0,0,1);
  pop.mut_lookup.insert(0.1);
  pop.mcounts.emplace_back(1);
  pop.mutations.emplace_back(0.2,0,0,1);
  pop.mut_lookup.insert(0.2);
  pop.mcounts.emplace_back(0);
  KTfwd::update_mutations(pop.mutations,pop.mut_lookup,pop.mcounts);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],1);
  BOOST_REQUIRE_EQUAL(pop.mcounts[1],0);
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),2);
  BOOST_REQUIRE_EQUAL(pop.fixations.size(),0);
  BOOST_REQUIRE_EQUAL(pop.fixation_times.size(),0);
}

BOOST_AUTO_TEST_CASE( test_recycle_all_fixations )
{
  unsigned N=1000;
  singlepop_t pop(N);
  //This is a neutral mutation
  pop.mutations.emplace_back(0.1,0,0,1);
  pop.mut_lookup.insert(0.1);
  pop.mcounts.emplace_back(2*N);
  KTfwd::update_mutations(pop.mutations,pop.mut_lookup,pop.mcounts,2*N);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],0);
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.fixations.size(),0);
  BOOST_REQUIRE_EQUAL(pop.fixation_times.size(),0);
}

BOOST_AUTO_TEST_CASE( test_recycle_all_fixations_and_move_to_fixations_vector )
{
  unsigned N=1000;
  singlepop_t pop(N);
  //This is a neutral mutation
  pop.mutations.emplace_back(0.1,0,0,1);
  pop.mut_lookup.insert(0.1);
  pop.mcounts.emplace_back(2*N);
  //use generation = 2
  KTfwd::update_mutations(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,2,2*N);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],0);
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.fixations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.fixation_times.size(),1);
  BOOST_REQUIRE_EQUAL(pop.fixation_times[0],2);
}


BOOST_AUTO_TEST_CASE( only_recycle_neutral_fixations )
{
  unsigned N=1000;
  singlepop_t pop(N);
  //This is a neutral mutation
  pop.mutations.emplace_back(0.1,0,0,1);
  pop.mut_lookup.insert(0.1);
  pop.mcounts.emplace_back(2*N);
  //This is NOT a neutral mutation
  pop.mutations.emplace_back(0.2,-1.0,0,1);
  pop.mut_lookup.insert(0.2);
  pop.mcounts.emplace_back(2*N);
  //use generation = 2
  KTfwd::update_mutations_n(pop.mutations,pop.fixations,pop.fixation_times,pop.mut_lookup,pop.mcounts,2,2*N);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],0);
  BOOST_REQUIRE_EQUAL(pop.mcounts[1],2*N);
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),2);
  BOOST_REQUIRE_EQUAL(pop.fixations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.fixation_times.size(),1);
  BOOST_REQUIRE_EQUAL(pop.fixation_times[0],2);
}
