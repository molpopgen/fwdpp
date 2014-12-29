#define BOOST_TEST_MODULE mutationTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <iostream>
#include <utility>
#include <list>
/*
  Testing move semantics with mutations is tricky.
  The typical data elements are trivial types (double, int, etc.),
  and thus a move may perform a copy if doing so is faster.

  Here, we put a vector into our mutation class and check that it actually moves.
  A real-world application of this would be if mutations keep track of their
  frequency trajectories over time.
*/
struct mut : public KTfwd::mutation_base
{
  std::vector<int> stuff;
  mut(const double & position, const unsigned & count, const bool & isneutral) :
    KTfwd::mutation_base(position,count,isneutral)
  {
  }
};

using gtype = KTfwd::gamete;

BOOST_AUTO_TEST_CASE( move_construct )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );

  mut m2(std::move(m1));

  BOOST_CHECK_EQUAL( m1.stuff.size(), 0 );
  BOOST_CHECK_EQUAL( m2.stuff.size(), 4 ); 
}


BOOST_AUTO_TEST_CASE( move_assign )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  
  mut m2 = std::move(m1);

  BOOST_CHECK_EQUAL( m1.stuff.size(), 0 );
  BOOST_CHECK_EQUAL( m2.stuff.size(), 4 ); 
}

//A "local" version of KTfwd::insert_at_end template.
//This version implements c+11 "perfect forwarding"
template<typename T, typename cT>
inline typename cT::iterator local_insert_at_end(  T && t,  cT * ct )
{
  return ct->insert(ct->end(),std::forward<T>(t));
}

struct local_insert_at_end2
{
  template<typename T, typename cT>
  inline typename cT::iterator operator()( T && t, cT * ct ) const
  {
    return ct->insert(ct->end(),std::forward<T>(t));
  }
};

//The next three tests check that perfect forwarding works

//Test for move into list
BOOST_AUTO_TEST_CASE( policy_test_1 )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  std::list<mut> mlist;

  local_insert_at_end( std::move(m1), &mlist );

  BOOST_CHECK_EQUAL( m1.stuff.size(), 0 );
}

//Test for copy into list
BOOST_AUTO_TEST_CASE( policy_test_2 )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  std::list<mut> mlist;

  local_insert_at_end( std::cref(m1), &mlist );

  BOOST_CHECK_EQUAL( m1.stuff.size(), 4 );
}

//Test for copy into list
BOOST_AUTO_TEST_CASE( policy_test_3 )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  std::list<mut> mlist;

  local_insert_at_end( std::ref(m1), &mlist );

  BOOST_CHECK_EQUAL( m1.stuff.size(), 4 );
}

//A "fake" mutation function that passes std::move to
//the policy
template<typename policy>
std::pair< std::list<mut>::iterator, bool >
faux_mutate( std::list<mut> * mlist,
	     const policy & mpolicy )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  auto itr = mpolicy(std::move(m1),mlist);
  return std::make_pair(itr,m1.stuff.empty());
}

//A "fake" mutation function that passes const reference
//the policy
template<typename policy>
std::pair< std::list<mut>::iterator, bool >
faux_mutate2( std::list<mut> * mlist,
	      const policy & mpolicy )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  auto itr = mpolicy(std::cref(m1),mlist);
  return std::make_pair(itr,m1.stuff.empty());
}

//Does passing stuff through a policy still get "perfectly forwarded"?
BOOST_AUTO_TEST_CASE( policy_test_4 )
{
  std::list<mut> mlist;
  auto rv = faux_mutate( &mlist,
			 std::bind( local_insert_at_end2(),std::placeholders::_1,std::placeholders::_2));
  BOOST_CHECK_EQUAL(rv.first->stuff.size(),4);
  BOOST_CHECK_EQUAL(rv.second,true);
}

BOOST_AUTO_TEST_CASE( policy_test_5 )
{
  std::list<mut> mlist;
  auto rv = faux_mutate2( &mlist,
			  std::bind( local_insert_at_end2(),std::placeholders::_1,std::placeholders::_2));
  BOOST_CHECK_EQUAL(rv.first->stuff.size(),4);
  BOOST_CHECK_EQUAL(rv.second,false);
}

//This is tricky: note what you have to do in the policy binding...
BOOST_AUTO_TEST_CASE( policy_test_6 )
{
  std::list<mut> mlist;
  auto rv = faux_mutate2( &mlist,
			  std::bind( local_insert_at_end<mut const &,std::list<mut> >,std::placeholders::_1,std::placeholders::_2));
  BOOST_CHECK_EQUAL(rv.first->stuff.size(),4);
  BOOST_CHECK_EQUAL(rv.second,false);
}

//Lambda expressions for policies
BOOST_AUTO_TEST_CASE( policy_test_7 )
{
  std::list<mut> mlist;
  auto rv = faux_mutate( &mlist,
			 []( mut && m, std::list<mut> * __mlist ) {
			   return __mlist->insert(__mlist->end(),std::move(m));
			 });
  BOOST_CHECK_EQUAL(rv.first->stuff.size(),4);
  BOOST_CHECK_EQUAL(rv.second,true);
}

BOOST_AUTO_TEST_CASE( policy_test_8 )
{
  std::list<mut> mlist;
  auto rv = faux_mutate( &mlist,
			 []( const mut & m, std::list<mut> * __mlist ) {
			   return __mlist->insert(__mlist->end(),m);
			 });
  BOOST_CHECK_EQUAL(rv.first->stuff.size(),4);
  BOOST_CHECK_EQUAL(rv.second,false);
}

BOOST_AUTO_TEST_CASE( policy_test_9 )
{
  std::list<mut> mlist;
  auto rv = faux_mutate2( &mlist,
			 []( const mut & m, std::list<mut> * __mlist ) {
			   return __mlist->insert(__mlist->end(),m);
			 });
  BOOST_CHECK_EQUAL(rv.first->stuff.size(),4);
  BOOST_CHECK_EQUAL(rv.second,false);
}


//Great, so what does the built-in policy do?

//When confronted with a move, it should copy but not move, as it takes const mtype &
BOOST_AUTO_TEST_CASE( fwdpp_policy_test_1 )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  std::list<mut> mlist;

  KTfwd::insert_at_end( std::move(m1), &mlist );

  BOOST_CHECK_EQUAL( m1.stuff.size(), 4 );
  BOOST_CHECK_EQUAL( mlist.empty(), false );
}

BOOST_AUTO_TEST_CASE( fwdpp_policy_test_2 )
{
  mut m1(0.123,1,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  std::list<mut> mlist;

  auto rv = faux_mutate( &mlist,
			 std::bind( KTfwd::insert_at_end<mut,std::list<mut> >,std::placeholders::_1,std::placeholders::_2));
  BOOST_CHECK_EQUAL(mlist.empty(),false);
  BOOST_CHECK_EQUAL(rv.first->stuff.size(),4);
  BOOST_CHECK_EQUAL(rv.second,false);
}
