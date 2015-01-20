//! \file mlocusCrossoverTest.cc \ingroup unit

#define BOOST_TEST_MODULE mlocusCrossoverTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <list>
#include <iostream>
#include <gsl/gsl_rng.h>

using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;

/*
  Initiate random number generation system -- 
  these tests will not be random, but these objects
  are req'd for function calls to fwdpp's internal 
  stuff 
*/
gsl_rng * r =  gsl_rng_alloc(gsl_rng_ranlxs2);

//Set up our basic containers
using glist = std::list<gtype>;
using gvector = std::vector< glist >;
using mutlist = std::list<mut>;
using diploid_t = std::vector< std::pair< glist::iterator, glist::iterator> >;

/*
  Set up the following config:
  g1,l1 = 0.5
  g1,l2 = 0.75,
  g2,l1 = 0.25,
  g2,l2 = 0.9

  This function is all setup. 
  Typically, these steps would be random outcomes
  of a model, but we need to manually do these
  steps ourselves for unit testing.  In other words,
  what happens in this function is NOT what is happening
  within fwdpp.  However, the unit test modules below
  use what this function does as input to fwdpp's internal
  functions.
*/
void setup1( gvector & gametes,
	     mutlist & mlist,
	     diploid_t & diploid )
{
  gametes = gvector( 2,glist(2,gtype(1)) );
  mlist = mutlist();

  //To set this up, let's add the mutations:
  auto m1 = mlist.insert(mlist.end(),mut(0.5,0.,1));
  auto m2 = mlist.insert(mlist.end(),mut(0.75,0.,1));
  auto m3 = mlist.insert(mlist.end(),mut(0.25,0.,1));
  auto m4 = mlist.insert(mlist.end(),mut(0.9,0.,1));

  //Add the mutations to the gametes
  auto gitr = gametes[0].begin();
  gitr->mutations.push_back(m1);
  ++gitr;
  gitr->mutations.push_back(m3);
  gitr = gametes[1].begin();
  gitr->mutations.push_back(m2);
  ++gitr;
  gitr->mutations.push_back(m4);

  //Now, make a diploid
  diploid = diploid_t(2);
  gitr=gametes[0].begin();
  diploid[0].first = gitr;
  ++gitr;
  diploid[0].second = gitr;
  gitr = gametes[1].begin();
  diploid[1].first = gitr;
  ++gitr;
  diploid[1].second = gitr;
}

BOOST_AUTO_TEST_CASE( two_locus_test_1 )
{
  gvector gametes;
  mutlist mlist;
  diploid_t diploid;
  setup1(gametes,mlist,diploid);

  //copy!
  auto dip2(diploid);
  //Make a pointer to the 2nd locus of dip2
  auto ptr2cdip = dip2.begin()+1;

  //There WILL be a crossover
  double r_bw_loci = 1;

  bool p1 = true, LO = false;
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
  ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
							  //No rec w/in loci
							  [](glist::iterator &a,glist::iterator & b) { return 0; },
							  //Rec. b/w loci returns an ODD number
							  [](gsl_rng * __r, const double & __d) { return 1; },
							  &r_bw_loci,1,
							  //the parental gamete types
							  diploid[1].first,diploid[1].second,
							  p1,LO);

  //Now, there has been a swap
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.9 );
}

BOOST_AUTO_TEST_CASE( two_locus_test_2 )
{
  //Redo it so that there is not a swap
  gvector gametes;
  mutlist mlist;
  diploid_t diploid;
  setup1(gametes,mlist,diploid);

  //copy!
  auto dip2(diploid);
  //Make a pointer to the 2nd locus of dip2
  auto ptr2cdip = dip2.begin()+1;

  //There WILL be a crossover
  double r_bw_loci = 1;
  
  bool p1 = true;
  bool LO = false;
  ptr2cdip = dip2.begin()+1;
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
  ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
  							  //No rec w/in loci
  							  [](glist::iterator &a,glist::iterator & b) { return 0; },
  							  //Rec. b/w loci returns an EVEN number
  							  [](gsl_rng * __r, const double & __d) { return 0; },
  							  &r_bw_loci,1,
  							  //the parental gamete types
  							  diploid[1].first,diploid[1].second,
  							  p1,LO);

  //Now, there has NOT been a swap
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
}

BOOST_AUTO_TEST_CASE( two_locus_test_3 )
{
  //Redo the first example, with a rec event in locus 2, and an odd no. recs b/w loci
  //THIS TEST WILL BE HARDER TO DO!!!!
  gvector gametes;
  mutlist mlist;
  diploid_t diploid;
  setup1(gametes,mlist,diploid);

  //There WILL be a crossover
  double r_bw_loci = 1;
  auto dip2 = diploid;
  bool p1 = true;
  bool LO = false;
  auto ptr2cdip = dip2.begin()+1;
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
  ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
  							  //No Rec. rate = 1
							  std::bind(KTfwd::genetics101(),std::placeholders::_1,std::placeholders::_2,&gametes[1],1.,r,
								    []() { return 0.8;} ),
  							  //Rec. b/w loci returns an ODD number
  							  [](gsl_rng * __r, const double & __d) { return 0; },
  							  &r_bw_loci,1,
  							  //the parental gamete types
  							  diploid[1].first,diploid[1].second,
  							  p1,LO);
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
}
