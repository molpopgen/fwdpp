/*!
  \file gametaTest.cc
  \ingroup unit
  \brief Tests construction and assigment to gametes via std::move
*/
#define BOOST_TEST_MODULE gameteTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <iostream>
#include <vector>

using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;

BOOST_AUTO_TEST_CASE( move_construct )
{
  //Neutral mutations at positions 0.1 and 0.9, resp.
  gtype g1(1),g2(1);
  std::vector<mut> mvector(1,mut(0.1,0.));

  KTfwd::fwdpp_internal::add_new_mutation(0,mvector,g1);
  mvector.emplace_back(0.9,0.);
  KTfwd::fwdpp_internal::add_new_mutation(1,mvector,g2);

  gtype g3( std::move(g2) );
  BOOST_CHECK_EQUAL( g1.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 0 );
  BOOST_CHECK_EQUAL( g3.mutations.size(), 1 );
}


BOOST_AUTO_TEST_CASE( move_assign )
{
  //Neutral mutations at positions 0.1 and 0.9, resp.
  gtype g1(1),g2(1);
  std::vector<mut> mvector(1,mut(0.1,0.));

  KTfwd::fwdpp_internal::add_new_mutation(0,mvector,g1);
  mvector.emplace_back(0.9,0.);
  KTfwd::fwdpp_internal::add_new_mutation(1,mvector,g2);


  gtype g3 = std::move(g2);
  BOOST_CHECK_EQUAL( g1.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g2.mutations.size(), 0 );
  BOOST_CHECK_EQUAL( g3.mutations.size(), 1 );
}
