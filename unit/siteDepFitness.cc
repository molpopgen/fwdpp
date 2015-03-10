/*
  Unit tests of KTfwd::site_dependent_fitness
 */

#define BOOST_TEST_MODULE siteDepFitness
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <iostream>
#include <cmath>
#include <list>

using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;
using mlist = std::list<mut>;
using glist = std::list<gtype>;

/*
  This test creates a situation
  where gamete1 has a selected mutation and gamete2 does not.
  If issue #8 were to cause fitness to be mis-calculated,
  then this test will fail.

  However, it does not.  Even with the bug, the remaining bit of the function
  gets the calculation right.  Yay!
 */
BOOST_AUTO_TEST_CASE( simple_multiplicative1 )
{
  gtype g1(1),g2(1);
  mlist mutations;

  //add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use the dominance...)
  auto mitr = mutations.insert( mutations.end(), mut(0.1,0.1,1u,0.5) );
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g1);
  BOOST_CHECK_EQUAL(g1.smutations.size(),1);

  glist g;
  auto gitr1 = g.insert(g.end(),g1);
  auto gitr2 = g.insert(g.end(),g2);

  double w = KTfwd::site_dependent_fitness()(gitr1,gitr2,
					     [&](double & fitness,const mlist::iterator & __mut)
					     {
					       fitness *= std::pow(1. + __mut->s,2.);
					     },
					     [](double & fitness,const mlist::iterator & __mut)
					     {
					       fitness *= (1. + __mut->s);
					     },
					     1.);
  BOOST_CHECK_EQUAL(w,1.1);
}

/*
  g2 has it, g1 does not
 */
BOOST_AUTO_TEST_CASE( simple_multiplicative2 )
{
  gtype g1(1),g2(1);
  mlist mutations;

  //add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use the dominance...)
  auto mitr = mutations.insert( mutations.end(), mut(0.1,0.1,1u,0.5) );
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g2);
  BOOST_CHECK_EQUAL(g1.smutations.size(),0);
  BOOST_CHECK_EQUAL(g2.smutations.size(),1);

  glist g;
  auto gitr1 = g.insert(g.end(),g1);
  auto gitr2 = g.insert(g.end(),g2);

  double w = KTfwd::site_dependent_fitness()(gitr1,gitr2,
					     [&](double & fitness,const mlist::iterator & __mut)
					     {
					       fitness *= std::pow(1. + __mut->s,2.);
					     },
					     [](double & fitness,const mlist::iterator & __mut)
					     {
					       fitness *= (1. + __mut->s);
					     },
					     1.);
  BOOST_CHECK_EQUAL(w,1.1);
}

/*
  Both have it
*/
BOOST_AUTO_TEST_CASE( simple_multiplicative3 )
{
  gtype g1(1),g2(1);
  mlist mutations;

  //add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use the dominance...)
  auto mitr = mutations.insert( mutations.end(), mut(0.1,0.1,1u,0.5) );
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g2);
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g1);
  BOOST_CHECK_EQUAL(g1.smutations.size(),1);
  BOOST_CHECK_EQUAL(g2.smutations.size(),1);

  glist g;
  auto gitr1 = g.insert(g.end(),g1);
  auto gitr2 = g.insert(g.end(),g2);

  double w = KTfwd::site_dependent_fitness()(gitr1,gitr2,
					     [&](double & fitness,const mlist::iterator & __mut)
					     {
					       fitness *= std::pow(1. + __mut->s,2.);
					     },
					     [](double & fitness,const mlist::iterator & __mut)
					     {
					       fitness *= (1. + __mut->s);
					     },
					     1.);
  BOOST_CHECK_EQUAL(w,1.1*1.1);
}
