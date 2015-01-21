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
  //Redo the first example, with a rec event in locus 2, and an even no. recs b/w loci
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
							  [&gametes]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,0.8);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[1],g1,g2);
							  },
  							  //Rec. b/w loci returns an EVEN number
  							  [](gsl_rng * __r, const double & __d) { return 0; },
  							  &r_bw_loci,1,
  							  //the parental gamete types
  							  diploid[1].first,diploid[1].second,
  							  p1,LO);
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 2 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[1]->pos,0.9 );
}

BOOST_AUTO_TEST_CASE( two_locus_test_4 )
{
  //Redo the first example, with a rec event in locus 2, and an even no. recs b/w loci
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
 							  [&gametes]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,0.8);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[1],g1,g2);
							  },
  							  //Rec. b/w loci returns an ODD number
  							  [](gsl_rng * __r, const double & __d) { return 1; },
  							  &r_bw_loci,1,
  							  //the parental gamete types
  							  diploid[1].first,diploid[1].second,
  							  p1,LO);
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 0 );
  BOOST_CHECK_EQUAL( ptr2cdip->second->mutations.size(), 1 );
}

//A redo of 3, but with two diploids as parents, and we do the xover for both chromos
BOOST_AUTO_TEST_CASE( two_locus_test_5 )
{
  //Redo the first example, with a rec event in locus 2, and an even no. recs b/w loci
  gvector gametes;
  mutlist mlist;
  diploid_t diploid;
  setup1(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //this is the other parent
  //There WILL be a crossover
  double r_bw_loci = 1;
  auto dip2 = diploid;
  bool p1 = true;
  bool LO = false;
  auto ptr2cdip = dip2.begin()+1;
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );

  auto pcopy(diploid);
  ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&gametes]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,0.8);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[1],g1,g2);
							  },
  							  //Rec. b/w loci returns an EVEN number
  							  [](gsl_rng * __r, const double & __d) { return 0; },
  							  &r_bw_loci,1,
  							  //the parental gamete types
  							  pcopy[1].first,pcopy[1].second,
  							  p1,LO);
  BOOST_CHECK_EQUAL( diploid[1].first->mutations.size(),1 );
  BOOST_CHECK_EQUAL( diploid[1].second->mutations.size(),1 );
  p1=false; //set it equal to false, for fun
  LO=false;
  pcopy = diploid;
  ptr2cdip->second = KTfwd::fwdpp_internal::multilocus_rec(r,
 							  [&gametes]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,0.8);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[1],g1,g2);
							  },
							   //Rec. b/w loci returns an EVEN number
							   [](gsl_rng * __r, const double & __d) { return 0; },
							   &r_bw_loci,1,
							   //the parental gamete types
							   pcopy[1].first,pcopy[1].second,
							   p1,LO);
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 2 );
  BOOST_CHECK_EQUAL( ptr2cdip->second->mutations.size(), 0 );
}


//OK, time now for a "real" test:
BOOST_AUTO_TEST_CASE( two_locus_test_6 )
{
  gvector gametes;
  mutlist mlist;
  diploid_t diploid;
  setup1(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //create the second parent

  //diploid, and diploid2 are the parents
  diploid_t offspring(2);

  bool p1g1 = true,LO = false,
    p2g1 = true,L1 = false;

  double r_bw_loci = 1.;
  //manipulae the first locus, i = 0.
  unsigned i = 0;
  auto ptr2cdip = offspring.begin()+i;
  /*
    There is a crossover in the first parent before the first mutation at position 0.1, which will swap the two mutations in the parent
    from p_gam_1 to p_gam_2, and there WILL be an xover b/w the two loci in parent 1
  */
  ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&gametes,&i]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,0.1);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[i],g1,g2);
							  },
							  //Rec. b/w loci returns an ODD number, but this is not relevant b/c i = 0
							  [](gsl_rng * __r, const double & __d) { return 1; },
							  &r_bw_loci,i,
							  //the parental gamete types
							  diploid[i].first,diploid[i].second,
							  p1g1,LO);

  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 1 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos, 0.25 );
  BOOST_CHECK_EQUAL( p1g1, true );
  BOOST_CHECK_EQUAL( LO, true );
  /*
    In parent 2, there are two x-overs: one before the two mutations, and one in b/w the loci
   */

  ptr2cdip->second = KTfwd::fwdpp_internal::multilocus_rec(r,
							   [&gametes,&i]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,0.1);
							    pos.push_back(0.3);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[i],g1,g2);
							   },
							   //Rec. b/w loci returns an ODD number, but this is not relevant b/c i = 0
							   [](gsl_rng * __r, const double & __d) { return 1; },
							   &r_bw_loci,i,
							   //the parental gamete types come from parent 2
							   diploid2[i].first,diploid2[i].second,
							   p2g1,L1);

  BOOST_CHECK_EQUAL( ptr2cdip->second->mutations.size(), 2 );
  BOOST_CHECK_EQUAL( ptr2cdip->second->mutations[0]->pos, 0.25 ); 
  BOOST_CHECK_EQUAL( ptr2cdip->second->mutations[1]->pos, 0.5 ); 
  BOOST_CHECK_EQUAL( p2g1, true );
  BOOST_CHECK_EQUAL( L1, false );

  //Now, the second locus
  ++i;
  ++ptr2cdip;

  /*
    There is a crossover in locus 2 b/w the two mutations, AND a crossover beteen locus 1 and 2, in parent 1
  */
  ptr2cdip->first =  KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&gametes,&i]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,0.8);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[i],g1,g2);
							  },
							  //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							  [](gsl_rng * __r, const double & __d) { return 1; },
							  &r_bw_loci,i,
							  //the parental gamete types
							  diploid[i].first,diploid[i].second,
							  p1g1,LO);
  
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 2 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos, 0.75 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[1]->pos, 0.9 );
  BOOST_CHECK_EQUAL( p1g1, true );
  BOOST_CHECK_EQUAL( LO, true );
  /*
    No crossover in locus 2, no crossover b/w locus 1 and 2, in parent 2
   */
  ptr2cdip->second =  KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&gametes,&i]( glist::iterator & g1, glist::iterator & g2 ) {
							    std::vector<double> pos(1,std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes[i],g1,g2);
							  },
							  //Rec. b/w loci returns an EVEN number, which will cause NO x-over b/w loci 1 and 2
							  [](gsl_rng * __r, const double & __d) { return 0; },
							  &r_bw_loci,i,
							  //the parental gamete types
							  diploid2[i].first,diploid2[i].second,
							  p2g1,L1);
}

/*
  Setup fxn for 3-locus scenario

  Set up the following config:
  g1,l1 = 0.5
  g1,l2 = 0.75,
  g2,l1 = 0.25,
  g2,l2 = 0.9
  g1,l3 = 1.25, 1.5
  g2,l3 = 1.1
*/
void setup2( gvector & gametes,
	     mutlist & mlist,
	     diploid_t & diploid )
{
  gametes = gvector( 3,glist(2,gtype(1)) );
  mlist = mutlist();

  //To set this up, let's add the mutations:
  auto m1 = mlist.insert(mlist.end(),mut(0.5,0.,1));
  auto m2 = mlist.insert(mlist.end(),mut(0.75,0.,1));
  auto m3 = mlist.insert(mlist.end(),mut(0.25,0.,1));
  auto m4 = mlist.insert(mlist.end(),mut(0.9,0.,1));
  auto m5 = mlist.insert(mlist.end(),mut(1.25,0.1,1));
  auto m6 = mlist.insert(mlist.end(),mut(1.5,0.1,1));
  auto m7 = mlist.insert(mlist.end(),mut(1.1,0.1,1));

  //put the mutations into gametes
  auto gitr = gametes[0].begin(); //gamete 1, locus 1
  gitr->mutations.push_back(m1);
  ++gitr;                         //gamete 2, locus 1
  gitr->mutations.push_back(m3);
  gitr = gametes[1].begin();      //gamete 1, locus 2
  gitr->mutations.push_back(m2);
  ++gitr;                         //gamete 2, locus 2
  gitr->mutations.push_back(m4);
  gitr = gametes[2].begin();      //gamete 1, locus 3
  gitr->mutations.push_back(m5);
  gitr->mutations.push_back(m6);
  ++gitr;                         //gamete 2, locus 3
  gitr->mutations.push_back(m7);

   //Now, make a diploid
  diploid = diploid_t(3);
  gitr=gametes[0].begin();
  diploid[0].first = gitr;
  ++gitr;
  diploid[0].second = gitr;
  gitr = gametes[1].begin();
  diploid[1].first = gitr;
  ++gitr;
  diploid[1].second = gitr;
  gitr = gametes[2].begin();
  diploid[2].first = gitr;
  ++gitr;
  diploid[2].second = gitr;
}

BOOST_AUTO_TEST_CASE( three_locus_test_1 )
{
  gvector gametes;
  mutlist mlist;
  diploid_t diploid;
  setup2(gametes,mlist,diploid);

  //This block makes sure that setup2 is working as far as gametes/mutations:
  for( auto & gg : gametes )
    {
      for (auto & g : gg)
	{
	  for( auto & m : g.mutations ) std::cerr << m->pos << ' ';
	  std::cerr << " | ";
	}
      std::cerr << '\n';
    }

  //And the diploid:
  for( auto & d : diploid )
    {
      for( auto & m : d.first->mutations ) std::cerr << m->pos << ' ';
      std::cerr << " | ";
      for( auto & m : d.second->mutations ) std::cerr << m->pos << ' ';
      std::cerr << '\n';
    }
}
