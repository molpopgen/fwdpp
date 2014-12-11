#define BOOST_TEST_MODULE mutateTest
#define BOOST_TEST_DYN_LINK 

#include <fwdpp/diploid.hh>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>


//trivial ways to play with the KTfwd::mutation type
using mut = KTfwd::mutation;

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
