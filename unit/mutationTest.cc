/*!
  \file mutationTest.cc 
  \ingroup unit
  \brief Tests that move semantics work with mutations
*/
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
  mut(const double & position,const bool & isneutral) :
    KTfwd::mutation_base(position,isneutral)
  {
  }
};

using gtype = KTfwd::gamete;

BOOST_AUTO_TEST_CASE( move_construct )
{
  mut m1(0.123,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );

  mut m2(std::move(m1));

  BOOST_CHECK_EQUAL( m1.stuff.size(), 0 );
  BOOST_CHECK_EQUAL( m2.stuff.size(), 4 ); 
}


BOOST_AUTO_TEST_CASE( move_assign )
{
  mut m1(0.123,1);
  m1.stuff = std::vector<int>( {2,3,4,5} );
  
  mut m2 = std::move(m1);

  BOOST_CHECK_EQUAL( m1.stuff.size(), 0 );
  BOOST_CHECK_EQUAL( m2.stuff.size(), 4 ); 
}

