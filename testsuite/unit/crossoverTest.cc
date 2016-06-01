/*!
  \defgroup unit Unit testing

  This documentation page is a list of links to the source code for the library's unit test programs.

  These sources are distributed with fwdpp in the unit/ directory of the source repository.

  They are implemented using boost's unit testing framework, and will not be compiled unless
  the configure script finds this library on your system.
 */

/*! \file crossoverTest.cc
  \ingroup unit
  \brief Tests of KTfwd::fwdpp_internal::recombine_gametes
*/

#include <config.h>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/internal/recombination_common.hpp>
#include <fwdpp/debug.hpp>
#include <boost/test/unit_test.hpp>

using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;

/*
  These tests model what goes on inside of the recombination routines.
*/

BOOST_AUTO_TEST_SUITE( test_crossing_over )

BOOST_FIXTURE_TEST_CASE( simple_test_1,standard_empty_single_deme_fixture )
{
  //Neutral mutations at positions 0.1 and 0.9, resp.
  gtype g1(1),g2(1);
  mutations.emplace_back(0.1,0.);

  g1.mutations.push_back(0);
  mutations.emplace_back(0.9,0.);
  g2.mutations.push_back(1);
  BOOST_CHECK_EQUAL( g1.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 1 );

  /*
    Let's make a single x-over at position 0.5
  */
  std::vector<double> rec_positions = {0.5, std::numeric_limits<double>::max()};

  //We need a "gamete pool"
  gametes = { g1, g2 };

  auto & MUT = mutations; //HACK so that lambda compiles...
  //Make sure that our iterators are cool
  BOOST_CHECK( std::find_if( g1.mutations.begin(),
			     g1.mutations.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.1;
			     } ) != g1.mutations.end() );
  BOOST_CHECK( std::find_if( g1.mutations.begin(),
			     g1.mutations.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.9;
			     } ) == g1.mutations.end() );
  BOOST_CHECK( std::find_if( g2.mutations.begin(),
			     g2.mutations.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.9;
			     } ) != g2.mutations.end() );
  BOOST_CHECK( std::find_if( g2.mutations.begin(),
			     g2.mutations.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.1;
			     } ) == g2.mutations.end() );

  BOOST_CHECK_EQUAL( g1.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 1 );

  //These are our containers for the recombinants

  gtype::mutation_container neutral,selected;
  KTfwd::fwdpp_internal::recombine_gametes( rec_positions,
					    0,1,
					    gametes,mutations,
					    neutral,selected);
  /*
    Now, neutral should have both mutations
  */
  BOOST_CHECK( std::find_if( neutral.begin(),
			     neutral.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.9;
			     } ) != neutral.end() );

  BOOST_CHECK( std::find_if( neutral.begin(),
			     neutral.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.1;
			     } ) != neutral.end() );

}

BOOST_FIXTURE_TEST_CASE( three_point_cross_1,standard_empty_single_deme_fixture )
{
  //g1: neutral muts at 0.1,0.5
  //g2: neutral muts at 0.9
  gtype g1(1),g2(1);
  mutations.emplace_back(0.1,0.);

  g1.mutations.push_back(0);
  mutations.emplace_back(0.5,0.);
  g1.mutations.push_back(1);
  mutations.emplace_back(0.9,0.);
  g2.mutations.push_back(2);
  BOOST_CHECK_EQUAL( g1.mutations.size(), 2 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 1 );

  std::vector<double> rec_positions = {0.25, 0.75, std::numeric_limits<double>::max()};

  //We need a "gamete pool"
  gametes = { g1, g2 };

  //Needed as of 0.3.3
  gtype::mutation_container neutral,selected;
  KTfwd::fwdpp_internal::recombine_gametes( rec_positions,
					    0,1,gametes,mutations,
					    neutral,selected );


  /*
    This was a double x-over.
    Neutral must therefore only contain 0.1
  */
  auto & MUT = mutations; //HACK so that lambda compiles...
  BOOST_CHECK( std::find_if( neutral.begin(),
			     neutral.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.1;
			     } ) != neutral.end());

  for( auto d : {0.5,0.9} )
    {
      BOOST_CHECK( std::find_if( neutral.begin(),
				 neutral.end(),
				 [d,&MUT]( std::size_t __m ) {
				   return MUT[__m].pos == d;
				 } ) == neutral.end() );

    }
}

/*
  A unit test to ensure that issue #26 can't happen again
*/
BOOST_FIXTURE_TEST_CASE( three_point_cross_2,standard_empty_single_deme_fixture )
{
  //g1: neutral muts at 0.1,0.5
  //g2: neutral muts at 0.9
  gtype g1(1),g2(1);
  mutations.push_back(mut(0.1,0.));
  g1.mutations.push_back(0);
  mutations.emplace_back(0.5,0.);
  g1.mutations.push_back(1);
  mutations.emplace_back(0.9,0.);
  g2.mutations.push_back(2);
  BOOST_CHECK_EQUAL( g1.mutations.size(), 2 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 1 );

  /*
    Here, we differ from above test in that an x-over position equals a mutation position.
    This tests the sanity of our upper_bound search in rec_gamete_updater
   */
  std::vector<double> rec_positions = {0.25, 0.5, std::numeric_limits<double>::max()};

  //We need a "gamete pool"
  gametes = { g1, g2 };

  KTfwd::fwdpp_internal::recombine_gametes( rec_positions,
					    0,1,gametes,mutations,
					    neutral,selected );


  auto & MUT = mutations; //HACK so that lambda compiles...
  /*
    This was a double x-over.
    Neutral must therefore only contain 0.1
  */
  BOOST_CHECK( std::find_if( neutral.begin(),
			     neutral.end(),
			     [&MUT]( std::size_t __m ) {
			       return MUT[__m].pos == 0.1;
			     } ) != neutral.end());

  for( auto d : {0.5,0.9} )
    {
      BOOST_CHECK( std::find_if( neutral.begin(),
				 neutral.end(),
				 [d,&MUT]( std::size_t __m ) {
				   return MUT[__m].pos == d;
				 } ) == neutral.end() );

    }
}

BOOST_AUTO_TEST_SUITE_END()
