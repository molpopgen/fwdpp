/*!
  \file test_sugar_change_neutral.cc

  \brief test KTfwd::change_neutral
*/

#define BOOST_TEST_MODULE test_sugar_add_mutation
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include <fwdpp/sugar/add_mutation.hpp>
#include <fwdpp/sugar/change_neutral.hpp>

using mutation_t = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_t>;
using singlepop_t = KTfwd::singlepop<mutation_t>;


BOOST_AUTO_TEST_CASE( test_change_neutral )
{
  singlepop_t pop(1000);
  KTfwd::add_mutation(pop,
		      //individuals where we want to place the mutation
		      {0,1,3,5,7,9},
		      /*
			gametes in each individual: 0 = .first, 1 = .second, 2 = .first and .second
			Thus, there should be 1+1+1+2+2+1=8 copies of the mutation in the population
		      */
		      {0,1,0,2,2,0},
		      //Parameters to pass on to create a new mutation
		      0.1,-0.1,1,0);
  BOOST_REQUIRE_EQUAL(KTfwd::check_sum(pop.gametes,2000),true);
  BOOST_REQUIRE_EQUAL(pop.gametes.size(),2); 
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],8);
  BOOST_REQUIRE_EQUAL(pop.mutations[0].neutral,false);
  //Change the mutation from selected to neutral
  KTfwd::change_neutral(pop,0);
  //Change it back
  KTfwd::change_neutral(pop,0);
}
