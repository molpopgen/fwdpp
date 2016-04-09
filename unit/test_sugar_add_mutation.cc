/*!
  \file test_sugar_add_mutation.cc

  \brief test KTfwd::add_mutation
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
#include <fwdpp/sugar/util.hpp>

using mutation_t = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_t>;
using singlepop_t = KTfwd::singlepop<mutation_t>;
using metapop_t = KTfwd::metapop<mutation_t>;
using multiloc_t = KTfwd::multiloc<mutation_t>;

BOOST_AUTO_TEST_CASE( test_add_mutation )
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
  BOOST_REQUIRE_EQUAL(pop.gametes.size(),9); //8 copies of mutation generated hap-hazardly will mean 9 total gametes in pop
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],8);
  BOOST_REQUIRE_EQUAL(pop.mutations[0].neutral,false);

  for( auto i : {0,3,9} ) //should have mutation on first gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].second].smutations.size(),0 );
    }
  for( auto i : {1} ) //should have mutation on second gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].first].smutations.size(),0 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].second].smutations.size(),1 );
    }
  for(auto i : {5,7}) // should have mutations on both gametes
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].second].smutations.size(),1 );
    }
}

BOOST_AUTO_TEST_CASE( test_add_mutation_from_object )
{
  singlepop_t pop(1000);
  mutation_t m(0.1,-0.1,1,0);
  KTfwd::add_mutation(pop,
		      //individuals where we want to place the mutation
		      {0,1,3,5,7,9},
		      /*
			gametes in each individual: 0 = .first, 1 = .second, 2 = .first and .second
			Thus, there should be 1+1+1+2+2+1=8 copies of the mutation in the population
		      */
		      {0,1,0,2,2,0},
		      //move it right into place:
		      std::move(m));
  BOOST_REQUIRE_EQUAL(KTfwd::check_sum(pop.gametes,2000),true);
  BOOST_REQUIRE_EQUAL(pop.gametes.size(),9); //8 copies of mutation generated hap-hazardly will mean 9 total gametes in pop
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],8);
  BOOST_REQUIRE_EQUAL(pop.mutations[0].neutral,false);

  for( auto i : {0,3,9} ) //should have mutation on first gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].second].smutations.size(),0 );
    }
  for( auto i : {1} ) //should have mutation on second gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].first].smutations.size(),0 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].second].smutations.size(),1 );
    }
  for(auto i : {5,7}) // should have mutations on both gametes
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i].second].smutations.size(),1 );
    }
}

BOOST_AUTO_TEST_CASE( test_add_mutation_metapop )
{
  metapop_t pop({1000,1000});
  const std::size_t DEME = 1; //we're gonna add the mutation into the second deme
  /*
    Note: for meta-pops, function takes vectors and vectors of vectors,
    so use via an initializer list has an extra bracket around it...
   */
  KTfwd::add_mutation(pop,
		      //deme index...
		      {DEME},
		      //individuals where we want to place the mutation
		      {{{0,1,3,5,7,9}}},
		      /*
			gametes in each individual: 0 = .first, 1 = .second, 2 = .first and .second
			Thus, there should be 1+1+1+2+2+1=8 copies of the mutation in the population
		      */
		      {{{0,1,0,2,2,0}}},
		      //For fun, pass in new mutation as a temporary
		      mutation_t(0.1,-0.1,1,0));
  BOOST_REQUIRE_EQUAL(KTfwd::check_sum(pop.gametes,4000),true);
  BOOST_REQUIRE_EQUAL(pop.gametes.size(),9); //8 copies of mutation generated hap-hazardly will mean 9 total gametes in pop
  BOOST_REQUIRE_EQUAL(pop.mutations.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts.size(),1);
  BOOST_REQUIRE_EQUAL(pop.mcounts[0],8);
  BOOST_REQUIRE_EQUAL(pop.mutations[0].neutral,false);

  for( auto i : {0,3,9} ) //should have mutation on first gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[DEME][i].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[DEME][i].second].smutations.size(),0 );
    }
  for( auto i : {1} ) //should have mutation on second gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[DEME][i].first].smutations.size(),0 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[DEME][i].second].smutations.size(),1 );
    }
  for(auto i : {5,7}) // should have mutations on both gametes
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[DEME][i].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[DEME][i].second].smutations.size(),1 );
    }
}

BOOST_AUTO_TEST_CASE( test_add_mutation_multiloc )
{
  multiloc_t pop(1000,10); //10-locus system
  const std::size_t LOCUS=7; //we'll add mutations into the 8-th locus
  KTfwd::add_mutation(pop,
		      //locus index...
		      LOCUS,
		      //individuals where we want to place the mutation
		      {0,1,3,5,7,9},
		      /*
			gametes in each individual: 0 = .first, 1 = .second, 2 = .first and .second
			Thus, there should be 1+1+1+2+2+1=8 copies of the mutation in the population
		      */
		      {0,1,0,2,2,0},
		      //params for new mutation
		      0.1,-0.1,1,0);
  BOOST_REQUIRE_EQUAL(KTfwd::check_sum(pop.gametes,2000),true);
  for( auto i : {0,3,9} ) //should have mutation on first gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i][LOCUS].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i][LOCUS].second].smutations.size(),0 );
    }
  for( auto i : {1} ) //should have mutation on second gamete only
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i][LOCUS].first].smutations.size(),0 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i][LOCUS].second].smutations.size(),1 );
    }
  for(auto i : {5,7}) // should have mutations on both gametes
    {
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i][LOCUS].first].smutations.size(),1 );
      BOOST_REQUIRE_EQUAL( pop.gametes[pop.diploids[i][LOCUS].second].smutations.size(),1 );
    }
}
