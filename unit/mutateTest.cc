#define BOOST_TEST_MODULE mutateTest
#define BOOST_TEST_DYN_LINK 

#include <fwdpp/diploid.hh>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <list>
//trivial ways to play with the KTfwd::mutation type
using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;

BOOST_AUTO_TEST_CASE( make_mutation1 )
{
  //Mutation at position 0.1, selection coefficient of 0, count in population of 1
  mut m(0.1,0.,1);

  BOOST_REQUIRE_EQUAL(m.pos,0.1);
  BOOST_REQUIRE_EQUAL(m.n,1);
  BOOST_REQUIRE_EQUAL(m.neutral,true);
}

BOOST_AUTO_TEST_CASE( make_mutation2 )
{
  //Mutation at position 0.1, selection coefficient of 0, count in population of 1
  //dominance of 0.25
  mut m(0.1,0.,1,0.25);

  BOOST_REQUIRE_EQUAL(m.pos,0.1);
  BOOST_REQUIRE_EQUAL(m.n,1);
  BOOST_REQUIRE_EQUAL(m.neutral,true);
  BOOST_REQUIRE_EQUAL(m.h,0.25);
}

BOOST_AUTO_TEST_CASE( copy_construct1 )
{
  mut m(0.1,0.,1);
  mut m2(m);

  BOOST_REQUIRE(m == m2);
}

BOOST_AUTO_TEST_CASE( assign1 )
{
  mut m(0.1,0.,1);
  mut m2 = m;

  BOOST_REQUIRE(m == m2);
}

//Tests the internal machinery for adding a mutation to a gamete
BOOST_AUTO_TEST_CASE( add_mutation1 )
{
  //create a gamete at frequency 1
  gtype g(1);

  //create a neutral mutation and add it to the "mutation pool"
  std::list<mut> mlist;
  auto mitr = mlist.insert(mlist.end(),mut(0.1,0.,1));

  KTfwd::fwdpp_internal::add_new_mutation(mitr,g);

  //So now, there is 1 neutral mutation and no selected mutations
  BOOST_CHECK_EQUAL( g.mutations.size(), 1 );
  BOOST_CHECK( g.smutations.empty() );

  //let's put in a selected mutation
  mitr = mlist.insert(mlist.end(),mut(0.1,-2.,1));
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g);

  BOOST_CHECK_EQUAL( g.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g.smutations.size(), 1 );
}
