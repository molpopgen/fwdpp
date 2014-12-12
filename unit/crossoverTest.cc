#define BOOST_TEST_MODULE crossoverTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <list>

using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;

/*
  These tests model what goes on inside of the recombination routines.
*/

BOOST_AUTO_TEST_CASE( simple_test_1 )
{
  //Neutral mutations at positions 0.1 and 0.9, resp.
  gtype g1(1),g2(1);
  std::list<mut> mlist;

  auto mitr = mlist.insert(mlist.end(),mut(0.1,0.,1));
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g1);
  mitr = mlist.insert(mlist.end(),mut(0.9,0.,1));
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g2);

  BOOST_CHECK_EQUAL( g1.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 1 );

  /*
    Let's make a single x-over at position 0.5
  */
  std::vector<double> rec_positions = {0.5, std::numeric_limits<double>::max()};

  //We need a "gamete pool"
  std::vector<gtype> gametes = { g1, g2 };

  //get pointers/iterators to our two existing gametes
  auto g1_itr = gametes.begin(),
    g2_itr = g1_itr+1;

  //Make sure that our iterators are cool
  BOOST_CHECK( std::find_if( g1_itr->mutations.begin(),
			     g1_itr->mutations.end(),
			     []( std::list<mut>::iterator __m ) {
			       return __m->pos == 0.1;
			     } ) != g1_itr->mutations.end() );
  BOOST_CHECK( std::find_if( g1_itr->mutations.begin(),
			     g1_itr->mutations.end(),
			     []( std::list<mut>::iterator __m ) {
			       return __m->pos == 0.0;
			     } ) == g1_itr->mutations.end() );
  BOOST_CHECK( std::find_if( g2_itr->mutations.begin(),
			     g2_itr->mutations.end(),
			     []( std::list<mut>::iterator __m ) {
			       return __m->pos == 0.9;
			     } ) != g2_itr->mutations.end() );
  BOOST_CHECK( std::find_if( g2_itr->mutations.begin(),
			     g2_itr->mutations.end(),
			     []( std::list<mut>::iterator __m ) {
			       return __m->pos == 0.1;
			     } ) == g2_itr->mutations.end() );

  BOOST_CHECK_EQUAL( g1.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 1 );

  //These are our containers for the recombinants
  gtype ng1(1),ng2(1);
  
  KTfwd::fwdpp_internal::recombine_gametes( rec_positions,
					    g1_itr, g2_itr,
					    ng1,ng2 );

  /*
    Now, ng1 must have both mutations,
    and ng2 must be empty
  */
  BOOST_CHECK( std::find_if( ng1.mutations.begin(),
			     ng1.mutations.end(),
			     []( std::list<mut>::iterator __m ) {
			       return __m->pos == 0.9;
			     } ) != ng1.mutations.end() );

  BOOST_CHECK( std::find_if( ng1.mutations.begin(),
			     ng1.mutations.end(),
			     []( std::list<mut>::iterator __m ) {
			       return __m->pos == 0.1;
			     } ) != ng1.mutations.end() );

  BOOST_CHECK( ng2.mutations.empty() );
}

BOOST_AUTO_TEST_CASE( three_point_cross_1 )
{
  //g1: neutral muts at 0.1,0.5
  //g2: neutral muts at 0.9
  gtype g1(1),g2(1);
  std::list<mut> mlist;

  auto mitr = mlist.insert(mlist.end(),mut(0.1,0.,1));
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g1);
  mitr = mlist.insert(mlist.end(),mut(0.5,0.,1));
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g1);
  mitr = mlist.insert(mlist.end(),mut(0.9,0.,1));
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g2);

  BOOST_CHECK_EQUAL( g1.mutations.size(), 2 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 1 );

  std::vector<double> rec_positions = {0.25, 0.75, std::numeric_limits<double>::max()};

  //We need a "gamete pool"
  std::vector<gtype> gametes = { g1, g2 };

  //get pointers/iterators to our two existing gametes
  auto g1_itr = gametes.begin(),
    g2_itr = g1_itr+1;

  gtype ng1(1),ng2(1);
  
  KTfwd::fwdpp_internal::recombine_gametes( rec_positions,
					    g1_itr, g2_itr,
					    ng1,ng2 );


  /*
    This was a double x-over.  ng1 must have a mutation at 0.1
    and ng2 must have the mutations at 0.5 and 0.9
  */
  BOOST_CHECK( std::find_if( ng1.mutations.begin(),
			     ng1.mutations.end(),
			     []( std::list<mut>::iterator __m ) {
			       return __m->pos == 0.1;
			     } ) != ng1.mutations.end() );

  for( auto d : {0.5,0.9} )
    {
      BOOST_CHECK( std::find_if( ng2.mutations.begin(),
				 ng2.mutations.end(),
				 [d]( std::list<mut>::iterator __m ) {
				   return __m->pos == d;
				 } ) != ng2.mutations.end() );

    }
}
