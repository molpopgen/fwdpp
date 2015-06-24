/*!
  \file mlocusCrossoverTest.cc
  \ingroup unit
  \brief Tests KTfwd::fwdpp_internal::multilocus_rec
*/

#define BOOST_TEST_MODULE mlocusCrossoverTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <iostream>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <list>
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
void setup1( glist & gametes,
	     mutlist & mlist,
	     diploid_t & diploid )
{
  gametes = glist(4,gtype(1));
  mlist = mutlist();

  //To set this up, let's add the mutations:
  auto m1 = mlist.insert(mlist.end(),mut(0.5,0.,1));
  auto m2 = mlist.insert(mlist.end(),mut(0.75,0.,1));
  auto m3 = mlist.insert(mlist.end(),mut(0.25,0.,1));
  auto m4 = mlist.insert(mlist.end(),mut(0.9,0.,1));

  //Add the mutations to the gametes
  auto gitr = gametes.begin();
  gitr->mutations.push_back(m1);
  ++gitr;
  gitr->mutations.push_back(m3);
  ++gitr;
  gitr->mutations.push_back(m2);
  ++gitr;
  gitr->mutations.push_back(m4);

  //Now, make a diploid
  diploid = diploid_t(2);
  gitr=gametes.begin();
  diploid[0].first = gitr;
  ++gitr;
  diploid[0].second = gitr;
  ++gitr;
  diploid[1].first = gitr;
  ++gitr;
  diploid[1].second = gitr;
}

// BOOST_AUTO_TEST_CASE( two_locus_test_1 )
// {
//   glist gametes;
//   mutlist mlist;
//   diploid_t diploid;
//   gtype::mutation_container neutral,selected; //req'd as of 0.3.3
//   setup1(gametes,mlist,diploid);

//   //copy!
//   auto dip2(diploid);
//   //Make a pointer to the 2nd locus of dip2
//   auto ptr2cdip = dip2.begin()+1;
//   //Check that we can copy
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );

//   //There WILL be a crossover
//   double r_bw_loci = 1;

//   bool p1 = true, LO = false, swapped = false;

//   auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);

//   ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
// 							  //No rec w/in loci
// 							  [](glist::iterator &a,glist::iterator & b,decltype(gamete_lookup) & ) { return 0; },
// 							  //Rec. b/w loci returns an ODD number
// 							  [](gsl_rng * __r, const double & __d) { return 1; },
// 							  &r_bw_loci,
// 							  1,
// 							  //the parental gamete types
// 							  diploid[1].first,
// 							  diploid[1].second,
// 							  //IDEA: The lookup to this locus' gametes.  Ought to streamline?
// 							  gamete_lookup,
// 							  p1,LO,swapped);

//   //Now, there has been a swap
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.9 );
// }

// BOOST_AUTO_TEST_CASE( two_locus_test_2 )
// {
//   //Redo it so that there is not a swap
//   glist gametes;
//   mutlist mlist;
//   diploid_t diploid;
//   setup1(gametes,mlist,diploid);

//   //copy!
//   auto dip2(diploid);
//   //Make a pointer to the 2nd locus of dip2
//   auto ptr2cdip = dip2.begin()+1;

//   //There WILL be a crossover
//   double r_bw_loci = 1;
  
//   bool p1 = true;
//   bool LO = false;
//   bool swapped = false;
//   ptr2cdip = dip2.begin()+1;
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
//   auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
//   ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
//   							  //No rec w/in loci
//   							  [](glist::iterator &a,glist::iterator & b,decltype(gamete_lookup) &) { return 0; },
//   							  //Rec. b/w loci returns an EVEN number
//   							  [](gsl_rng * __r, const double & __d) { return 0; },
//   							  &r_bw_loci,1,
//   							  //the parental gamete types
//   							  diploid[1].first,diploid[1].second,
// 							  gamete_lookup,
//   							  p1,LO,swapped);

//   //Now, there has NOT been a swap
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
// }

// BOOST_AUTO_TEST_CASE( two_locus_test_3 )
// {
//   //Redo the first example, with a rec event in locus 2, and an even no. recs b/w loci
//   glist gametes;
//   mutlist mlist;
//   diploid_t diploid;
//   setup1(gametes,mlist,diploid);

//   //There WILL be a crossover
//   double r_bw_loci = 1;
//   auto dip2 = diploid;
//   bool p1 = true;
//   bool LO = false;
//   bool swapped = false;
//   auto ptr2cdip = dip2.begin()+1;
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
//   BOOST_CHECK_EQUAL( ptr2cdip->second->mutations[0]->pos,0.9 );
//   auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
//   gtype::mutation_container neutral,selected; //req'd as of 0.3.3
//   ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
// 							  [&]( glist::iterator & g1, glist::iterator & g2,decltype(gamete_lookup) & ) {
// 							    std::vector<double> pos(1,0.8);
// 							    pos.push_back(std::numeric_limits<double>::max());
// 							    //Make use of overload that takes fixed number of positions instead of genetic map policy
// 							    return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
// 							  },
//   							  //Rec. b/w loci returns an EVEN number
//   							  [](gsl_rng * __r, const double & __d) { return 0; },
//   							  &r_bw_loci,1,
//   							  //the parental gamete types
//   							  diploid[1].first,diploid[1].second,
// 							  gamete_lookup,
//   							  p1,LO,swapped);
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 2 );
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[1]->pos,0.9 );
// }

// BOOST_AUTO_TEST_CASE( two_locus_test_4 )
// {
//   //Redo the first example, with a rec event in locus 2, and an even no. recs b/w loci
//   glist gametes;
//   mutlist mlist;
//   diploid_t diploid;
//   setup1(gametes,mlist,diploid);

//   //There WILL be a crossover
//   double r_bw_loci = 1;
//   auto dip2 = diploid;
//   bool p1 = true;
//   bool LO = false, swapped=false;
//   auto ptr2cdip = dip2.begin()+1;
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );
//   BOOST_CHECK_EQUAL( ptr2cdip->second->mutations[0]->pos,0.9 );
//   auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
//   gtype::mutation_container neutral,selected; //req'd as of 0.3.3
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 1 );
//   BOOST_CHECK_EQUAL( ptr2cdip->second->mutations.size(), 1 );
//   ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
//  							  [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & ) {
// 							    std::vector<double> pos(1,0.8);
// 							    pos.push_back(std::numeric_limits<double>::max());
// 							    //Make use of overload that takes fixed number of positions instead of genetic map policy
// 							    return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
// 							  },
//   							  //Rec. b/w loci returns an ODD number
//   							  [](gsl_rng * __r, const double & __d) { return 1; },
//   							  &r_bw_loci,1,
//   							  //the parental gamete types
//   							  diploid[1].first,diploid[1].second,
// 							  gamete_lookup,
//   							  p1,LO,swapped);
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 0 );
//   BOOST_CHECK_EQUAL( ptr2cdip->second->mutations.size(), 1 );
//   BOOST_CHECK_EQUAL( ptr2cdip->second->mutations[0]->pos,0.9 );
// }

// // //A redo of 3, but with two diploids as parents, and we do the xover for both chromos
// BOOST_AUTO_TEST_CASE( two_locus_test_5 )
// {
//   //Redo the first example, with a rec event in locus 2, and an even no. recs b/w loci
//   glist gametes;
//   mutlist mlist;
//   diploid_t diploid;
//   setup1(gametes,mlist,diploid);
//   diploid_t diploid2(diploid); //this is the other parent
//   //There WILL be a crossover
//   double r_bw_loci = 1;
//   auto dip2 = diploid;
//   bool p1 = true;
//   bool LO = false,swapped1=false,swapped2=false;
//   auto ptr2cdip = dip2.begin()+1;
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos,0.75 );

//   auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
//   gtype::mutation_container neutral,selected; //req'd as of 0.3.3

//   auto pcopy(diploid);
//   ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
// 							  [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & ) {
// 							    std::vector<double> pos(1,0.8);
// 							    pos.push_back(std::numeric_limits<double>::max());
// 							    //Make use of overload that takes fixed number of positions instead of genetic map policy
// 							    return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
// 							  },
//   							  //Rec. b/w loci returns an EVEN number
//   							  [](gsl_rng * __r, const double & __d) { return 0; },
//   							  &r_bw_loci,1,
//   							  //the parental gamete types
//   							  pcopy[1].first,pcopy[1].second,
// 							  gamete_lookup,
//   							  p1,LO,swapped1);
//   BOOST_CHECK_EQUAL( diploid[1].first->mutations.size(),1 );
//   BOOST_CHECK_EQUAL( diploid[1].second->mutations.size(),1 );
//   p1=false; //set it equal to false, for fun
//   LO=false;
//   pcopy = diploid;
//   for(auto mitr = pcopy[1].first->mutations.begin() ; mitr !=pcopy[1].first->mutations.end() ; ++mitr )
//     std::cerr << (*mitr)->pos << ' ';
//   std::cerr << '\n';

//   for(auto mitr = pcopy[1].second->mutations.begin() ; mitr !=pcopy[1].second->mutations.end() ; ++mitr )
//     std::cerr << (*mitr)->pos << ' ';
//   std::cerr << '\n';
//   ptr2cdip->second = KTfwd::fwdpp_internal::multilocus_rec(r,
// 							   [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & ) {
// 							    std::vector<double> pos(1,0.8);
// 							    pos.push_back(std::numeric_limits<double>::max());
// 							    //Make use of overload that takes fixed number of positions instead of genetic map policy
// 							    return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
// 							  },
// 							   //Rec. b/w loci returns an EVEN number
// 							   [](gsl_rng * __r, const double & __d) { return 0; },
// 							   &r_bw_loci,1,
// 							   //the parental gamete types
// 							   pcopy[1].first,pcopy[1].second,
// 							   gamete_lookup,
// 							   p1,LO,false);
//   for(auto mitr = pcopy[1].first->mutations.begin() ; mitr !=pcopy[1].first->mutations.end() ; ++mitr )
//     std::cerr << (*mitr)->pos << ' ';
//   std::cerr << '\n';

//   for(auto mitr = pcopy[1].second->mutations.begin() ; mitr !=pcopy[1].second->mutations.end() ; ++mitr )
//     std::cerr << (*mitr)->pos << ' ';
//   std::cerr << '\n';
//   BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 2 );
//   BOOST_CHECK_EQUAL( ptr2cdip->second->mutations.size(), 0 );
// }


// //OK, time now for a "real" test:
BOOST_AUTO_TEST_CASE( two_locus_test_6 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;
  setup1(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //create the second parent

  //diploid, and diploid2 are the parents
  diploid_t offspring(2);

  bool p1g1 = true,LO = false,
    p2g1 = true,L1 = false,swapped1=false,swapped2=false;

  double r_bw_loci = 1.;
  //manipulae the first locus, i = 0.
  unsigned i = 0;
  auto ptr2cdip = offspring.begin()+i;
  /*
    There is a crossover in the first parent before the first mutation at position 0.1, which will swap the two mutations in the parent
    from p_gam_1 to p_gam_2, and there WILL be an xover b/w the two loci in parent 1
  */
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  ptr2cdip->first = KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & ) {
							    std::vector<double> pos(1,0.1);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
							  },
							  //Rec. b/w loci returns an ODD number, but this is not relevant b/c i = 0
							  [](gsl_rng * __r, const double & __d) { return 1; },
							  &r_bw_loci,i,
							  //the parental gamete types
							  diploid[i].first,diploid[i].second,
							  gamete_lookup,
							  p1g1,LO,swapped1);

  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 1 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos, 0.25 );
  BOOST_CHECK_EQUAL( p1g1, true );
  BOOST_CHECK_EQUAL( LO, true );
  /*
    In parent 2, there are two x-overs: one before the two mutations, and one in b/w the loci
   */

  ptr2cdip->second = KTfwd::fwdpp_internal::multilocus_rec(r,
							   [&]( glist::iterator & g1, glist::iterator & g2,decltype(gamete_lookup) & ) {
							    std::vector<double> pos(1,0.1);
							    pos.push_back(0.3);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
							   },
							   //Rec. b/w loci returns an ODD number, but this is not relevant b/c i = 0
							   [](gsl_rng * __r, const double & __d) { return 1; },
							   &r_bw_loci,i,
							   //the parental gamete types come from parent 2
							   diploid2[i].first,diploid2[i].second,
							   gamete_lookup,
							   p2g1,L1,swapped2);

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
							   [&]( glist::iterator & g1, glist::iterator & g2,decltype(gamete_lookup) & ) {
							    std::vector<double> pos(1,0.8);
							    pos.push_back(std::numeric_limits<double>::max());
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
							   },
							   //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							   [](gsl_rng * __r, const double & __d) { return 1; },
							   &r_bw_loci,i,
							   //the parental gamete types
							   diploid[i].first,diploid[i].second,
							   gamete_lookup,
							   p1g1,LO,swapped1);
  
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations.size(), 2 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[0]->pos, 0.75 );
  BOOST_CHECK_EQUAL( ptr2cdip->first->mutations[1]->pos, 0.9 );
  BOOST_CHECK_EQUAL( p1g1, true );
  BOOST_CHECK_EQUAL( LO, true );
  /*
    No crossover in locus 2, no crossover b/w locus 1 and 2, in parent 2
   */
  ptr2cdip->second =  KTfwd::fwdpp_internal::multilocus_rec(r,
							    [&]( glist::iterator & g1, glist::iterator & g2,decltype(gamete_lookup) & ) {
							      std::vector<double> pos(1,std::numeric_limits<double>::max());
							      //Make use of overload that takes fixed number of positions instead of genetic map policy
							      return KTfwd::recombine_gametes(pos,&gametes,g1,g2,gamete_lookup,neutral,selected);
							    },
							    //Rec. b/w loci returns an EVEN number, which will cause NO x-over b/w loci 1 and 2
							    [](gsl_rng * __r, const double & __d) { return 0; },
							    &r_bw_loci,i,
							    //the parental gamete types
							    diploid2[i].first,diploid2[i].second,
							    gamete_lookup,
							    p2g1,L1,swapped2);
}

// /*
//   Setup fxn for 3-locus scenario

//   Set up the following config:
//   g1,l1 = 0.5
//   g1,l2 = 0.75,
//   g2,l1 = 0.25,
//   g2,l2 = 0.9
//   g1,l3 = 1.25, 1.5
//   g2,l3 = 1.1
// */
void setup2( glist & gametes,
	     mutlist & mlist,
	     diploid_t & diploid )
{
  gametes = glist(6,gtype(1));
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
  auto gitr = gametes.begin(); //gamete 1, locus 1
  gitr->mutations.push_back(m1);
  ++gitr;                         //gamete 2, locus 1
  gitr->mutations.push_back(m3);
  ++gitr;
  gitr->mutations.push_back(m2);
  ++gitr;                         //gamete 2, locus 2
  gitr->mutations.push_back(m4);
  ++gitr;
  gitr->mutations.push_back(m5);
  gitr->mutations.push_back(m6);
  ++gitr;                         //gamete 2, locus 3
  gitr->mutations.push_back(m7);

   //Now, make a diploid
  diploid = diploid_t(3);
  gitr=gametes.begin();
  diploid[0].first = gitr;
  ++gitr;
  diploid[0].second = gitr;
  ++gitr;
  diploid[1].first = gitr;
  ++gitr;
  diploid[1].second = gitr;
  ++gitr;
  diploid[2].first = gitr;
  ++gitr;
  diploid[2].second = gitr;
}

BOOST_AUTO_TEST_CASE( three_locus_test_1 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;            //parent 1
  setup2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2
  //This block makes sure that setup2 is working as far as gametes/mutations:
  // std::cerr << "gametes:\n";
  // for( auto & g : gametes )
  //   {
  //     for( auto & m : g.mutations ) std::cerr << m->pos << ' ';
  //     std::cerr << " | ";
  //     std::cerr << '\n';
  //   }

  // //And the diploid:
  // std::cerr << "diploids:\n";
  // for( auto & d : diploid )
  //   {
  //     for( auto & m : d.first->mutations ) std::cerr << m->pos << ' ';
  //     std::cerr << " | ";
  //     for( auto & m : d.second->mutations ) std::cerr << m->pos << ' ';
  //     std::cerr << '\n';
  //   }


  /*
    Make vectors of events regarding what happens withn and between each
    locus, for each parent.
  */
  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector<std::vector<double> > { std::vector<double>{ 0.3,MVAL },
						  std::vector<double>{ 0.55,0.8,MVAL },
						  std::vector<double>{1.3,MVAL}
  };
  auto rec2 = std::vector< std::vector<double> > { std::vector<double>{ 0.45, MVAL },
						   std::vector<double>{ MVAL },
						   std::vector<double>{ MVAL }
  };
  //recombinations b/w loci?
  std::vector<unsigned> bw1 = { 1,0 },bw2 = {0,1};
  std::vector<double> r_bw_loci = {1.,1.,1.};
  diploid_t offspring(3); 
  bool p1g1 = false,p2g1=true,LO1=false,LO2=false,s1=false,s2=false;
  auto ptr = offspring.begin();
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  for( unsigned i = 0 ; i < offspring.size() ; ++i,++ptr )
    {
      ptr->first =  KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							    //if ( rec1[i].empty() ) return;  
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(rec1[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							  },
							  //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							  [&bw1,&i](gsl_rng * __r, const double & __d) { return bw1[i-1]; },
							  &r_bw_loci[0],i,
							  //the parental gamete types
							  diploid[i].first,diploid[i].second,
							  gamete_lookup,
							  p1g1,LO1,s1);
      if ( s1 ) {
	for( auto itr = diploid.begin() + i + 1 ; itr < diploid.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      if(i==0)
	{
	  BOOST_REQUIRE_EQUAL(p1g1,false);
	  BOOST_REQUIRE_EQUAL(LO1,true);
	}
      neutral.clear();
      std::cerr << "before: ";
      for(auto mitr = diploid2[i].first->mutations.begin();mitr!=diploid2[i].first->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
      for(auto mitr = diploid2[i].second->mutations.begin();mitr!=diploid2[i].second->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
      ptr->second =  KTfwd::fwdpp_internal::multilocus_rec(r,
							   [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							     //if ( rec1[i].empty() ) return;  
							     //Make use of overload that takes fixed number of positions instead of genetic map policy
							     return KTfwd::recombine_gametes(rec2[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							   },
							   //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							   [&bw2,&i](gsl_rng * __r, const double & __d) { return bw2[i-1]; },
							   &r_bw_loci[0],i,
							   //the parental gamete types
							   diploid2[i].first,diploid2[i].second,
							   gamete_lookup,
							   p2g1,LO2,s2);
      if ( s2 ) {
	for( auto itr = diploid2.begin() + i + 1 ; itr < diploid2.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      std::cerr << "after: ";
      for(auto mitr =ptr->second->mutations.begin();mitr!=ptr->second->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
    }
  //Properties of the variables, etc., related to what happend with the offspring's FIRST gamete:
  BOOST_CHECK_EQUAL( LO1, true );
  BOOST_CHECK_EQUAL(offspring[0].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[1].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations[0]->pos,1.25);
  
  BOOST_CHECK_EQUAL(offspring[0].second->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations[0]->pos,0.75);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations[0]->pos,1.1);
}

BOOST_AUTO_TEST_CASE( three_locus_test_2 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;            //parent 1
  setup2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2


  /*
    Make vectors of events regarding what happens withn and between each
    locus, for each parent.
  */
  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector<std::vector<double> > { std::vector<double>{ 0.3,MVAL },
						  std::vector<double>{ 0.55,0.8,MVAL },
						  std::vector<double>{1.3,MVAL}
  };
  auto rec2 = std::vector< std::vector<double> > { std::vector<double>{ 0.45, MVAL },
						   std::vector<double>{ MVAL },
						   std::vector<double>{ MVAL }
  };
  //recombinations b/w loci?
  std::vector<unsigned> bw1 = { 1,0 }, bw2={0,1};
  std::vector<double> r_bw_loci = {1.,1.,1.};
  diploid_t offspring(3); 
  bool p1g1 = false,p2g1=true,LO1=false,LO2=false,s1=false,s2=false;
  auto ptr = offspring.begin();
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  for( unsigned i = 0 ; i < offspring.size() ; ++i,++ptr )
    {
      ptr->first =  KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							    //if ( rec1[i].empty() ) return;  
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(rec1[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							  },
							  //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							  [&bw1,&i](gsl_rng * __r, const double & __d) { return bw1[i-1]; },
							  &r_bw_loci[0],i,
							  //the parental gamete types
							  diploid[i].first,diploid[i].second,
							  gamete_lookup,
							  p1g1,LO1,s1);
      if ( s1 ) {
	for( auto itr = diploid.begin() + i + 1 ; itr < diploid.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      if(i==0)
	{
	  BOOST_REQUIRE_EQUAL(p1g1,false);
	  BOOST_REQUIRE_EQUAL(LO1,true);
	}
      neutral.clear();
      /*
      std::cerr << "before: ";
      for(auto mitr = diploid2[i].first->mutations.begin();mitr!=diploid2[i].first->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
      for(auto mitr = diploid2[i].second->mutations.begin();mitr!=diploid2[i].second->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
      */
      ptr->second =  KTfwd::fwdpp_internal::multilocus_rec(r,
							   [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							     //if ( rec1[i].empty() ) return;  
							     //Make use of overload that takes fixed number of positions instead of genetic map policy
							     return KTfwd::recombine_gametes(rec2[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							   },
							   //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							   [&bw2,&i](gsl_rng * __r, const double & __d) { return bw2[i-1]; },
							   &r_bw_loci[0],i,
							   //the parental gamete types
							   diploid2[i].first,diploid2[i].second,
							   gamete_lookup,
							   p2g1,LO2,s2);
      if ( s2 ) {
	for( auto itr = diploid2.begin() + i + 1 ; itr < diploid2.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      /*
      std::cerr << "after: ";
      for(auto mitr =ptr->second->mutations.begin();mitr!=ptr->second->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
      */
    }
  //Properties of the variables, etc., related to what happend with the offspring's FIRST gamete:
  BOOST_CHECK_EQUAL( LO1, true );
  BOOST_CHECK_EQUAL(offspring[0].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[1].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations[0]->pos,1.25);
  
  BOOST_CHECK_EQUAL(offspring[0].second->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations[0]->pos,0.75);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations[0]->pos,1.1);
}

BOOST_AUTO_TEST_CASE( three_locus_test_3 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;            //parent 1
  setup2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2


  /*
    Make vectors of events regarding what happens withn and between each
    locus, for each parent.
  */
  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector<std::vector<double> > { std::vector<double>{ 0.3,MVAL },
						  std::vector<double>{ 0.55,0.8,MVAL },
						  std::vector<double>{1.3,MVAL}
  };
  auto rec2 = std::vector< std::vector<double> > { std::vector<double>{ 0.1,0.45, MVAL },
						   std::vector<double>{ MVAL },
						   std::vector<double>{ MVAL }
  };
  //recombinations b/w all loci this time...
  std::vector<unsigned> bw1 = { 1,0 }, bw2={1,1};
  std::vector<double> r_bw_loci = {1.,1.,1.};
  diploid_t offspring(3); 
  bool p1g1 = false,p2g1=true,LO1=false,LO2=false,s1=false,s2=false;
  auto ptr = offspring.begin();
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  for( unsigned i = 0 ; i < offspring.size() ; ++i,++ptr )
    {
      ptr->first =  KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							    //if ( rec1[i].empty() ) return;  
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(rec1[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							  },
							  //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							  [&bw1,&i](gsl_rng * __r, const double & __d) { return bw1[i-1]; },
							  &r_bw_loci[0],i,
							  //the parental gamete types
							  diploid[i].first,diploid[i].second,
							  gamete_lookup,
							  p1g1,LO1,s1);
      if ( s1 ) {
	for( auto itr = diploid.begin() + i + 1 ; itr < diploid.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      if(i==0)
	{
	  BOOST_REQUIRE_EQUAL(p1g1,false);
	  BOOST_REQUIRE_EQUAL(LO1,true);
	  //BOOST_REQUIRE_EQUAL(s1,false);
	}
      neutral.clear();
      // std::cerr << "before: ";
      // for(auto mitr = diploid2[i].first->mutations.begin();mitr!=diploid2[i].first->mutations.end();++mitr)
      // 	std::cerr << (*mitr)->pos << ' ';
      // std::cerr << '\n';
      // for(auto mitr = diploid2[i].second->mutations.begin();mitr!=diploid2[i].second->mutations.end();++mitr)
	//	std::cerr << (*mitr)->pos << ' ';
      //std::cerr << '\n';
      ptr->second =  KTfwd::fwdpp_internal::multilocus_rec(r,
							   [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							     //if ( rec1[i].empty() ) return;  
							     //Make use of overload that takes fixed number of positions instead of genetic map policy
							     return KTfwd::recombine_gametes(rec2[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							   },
							   //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							   [&bw2,&i](gsl_rng * __r, const double & __d) { return bw2[i-1]; },
							   &r_bw_loci[0],i,
							   //the parental gamete types
							   diploid2[i].first,diploid2[i].second,
							   gamete_lookup,
							   p2g1,LO2,s2);
      //std::cerr << "S2 = " << s2 << '\n';
      if ( s2 ) {
	for( auto itr = diploid2.begin() + i + 1 ; itr < diploid2.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      // std::cerr << "after: ";
      // for(auto mitr =ptr->second->mutations.begin();mitr!=ptr->second->mutations.end();++mitr)
      // 	std::cerr << (*mitr)->pos << ' ';
      // std::cerr << '\n';
    }
  BOOST_CHECK_EQUAL( LO1, true );
  BOOST_CHECK_EQUAL(offspring[0].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[1].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations[0]->pos,1.25);
  
  BOOST_CHECK_EQUAL(offspring[0].second->mutations.size(),2);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations[0]->pos,0.9);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations.size(),2);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations[0]->pos,1.25);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations[1]->pos,1.5);
}

BOOST_AUTO_TEST_CASE( three_locus_test_4 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;            //parent 1
  setup2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2


  /*
    Make vectors of events regarding what happens withn and between each
    locus, for each parent.
  */
  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector<std::vector<double> > { std::vector<double>{ 0.3,MVAL },
						  std::vector<double>{ 0.55,0.8,MVAL },
						  std::vector<double>{1.3,MVAL}
  };
  //This one will be the doozy
  auto rec2 = std::vector< std::vector<double> > { std::vector<double>{ 0.1,0.45, 0.51, MVAL },
						   std::vector<double>{ 0.6,0.7,0.99, MVAL },
						   std::vector<double>{ 1,1.2,1.3,1.55,MVAL }
  };
  //recombinations b/w all loci this time...
  std::vector<unsigned> bw1 = { 1,0 }, bw2={1,1};
  std::vector<double> r_bw_loci = {1.,1.,1.};
  diploid_t offspring(3); 
  bool p1g1 = false,p2g1=true,LO1=false,LO2=false,s1=false,s2=false;
  auto ptr = offspring.begin();
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  std::cerr << "test 4:\n";
  for( unsigned i = 0 ; i < offspring.size() ; ++i,++ptr )
    {
      ptr->first =  KTfwd::fwdpp_internal::multilocus_rec(r,
							  [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							    //if ( rec1[i].empty() ) return;  
							    //Make use of overload that takes fixed number of positions instead of genetic map policy
							    return KTfwd::recombine_gametes(rec1[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							  },
							  //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							  [&bw1,&i](gsl_rng * __r, const double & __d) { return bw1[i-1]; },
							  &r_bw_loci[0],i,
							  //the parental gamete types
							  diploid[i].first,diploid[i].second,
							  gamete_lookup,
							  p1g1,LO1,s1);
      if ( s1 ) {
	for( auto itr = diploid.begin() + i + 1 ; itr < diploid.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      if(i==0)
	{
	  BOOST_REQUIRE_EQUAL(p1g1,false);
	  BOOST_REQUIRE_EQUAL(LO1,true);
	  //BOOST_REQUIRE_EQUAL(s1,false);
	}
      neutral.clear();
      std::cerr << "before: ";
      for(auto mitr = diploid2[i].first->mutations.begin();mitr!=diploid2[i].first->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
      for(auto mitr = diploid2[i].second->mutations.begin();mitr!=diploid2[i].second->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
      ptr->second =  KTfwd::fwdpp_internal::multilocus_rec(r,
							   [&]( glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) &  ) {
							     //if ( rec1[i].empty() ) return;  
							     //Make use of overload that takes fixed number of positions instead of genetic map policy
							     return KTfwd::recombine_gametes(rec2[i],&gametes,g1,g2,gamete_lookup,neutral,selected);
							   },
							   //Rec. b/w loci returns an ODD number, which will cause an x-over b/w loci 1 and 2
							   [&bw2,&i](gsl_rng * __r, const double & __d) { return bw2[i-1]; },
							   &r_bw_loci[0],i,
							   //the parental gamete types
							   diploid2[i].first,diploid2[i].second,
							   gamete_lookup,
							   p2g1,LO2,s2);
      std::cerr << "S2 = " << s2 << '\n';
      if ( s2 ) {
	for( auto itr = diploid2.begin() + i + 1 ; itr < diploid2.end() ; ++itr )
	  {
	    std::swap( itr->first, itr->second );
	  }
      }
      std::cerr << "after: ";
      for(auto mitr =ptr->second->mutations.begin();mitr!=ptr->second->mutations.end();++mitr)
	std::cerr << (*mitr)->pos << ' ';
      std::cerr << '\n';
    }
  BOOST_CHECK_EQUAL( LO1, true );
  BOOST_CHECK_EQUAL(offspring[0].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[1].first->mutations.size(),0);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[2].first->mutations[0]->pos,1.25);
  
  BOOST_CHECK_EQUAL(offspring[0].second->mutations.size(),2);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations.size(),1);
  BOOST_CHECK_EQUAL(offspring[1].second->mutations[0]->pos,0.75);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations.size(),2);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations[0]->pos,1.1);
  BOOST_CHECK_EQUAL(offspring[2].second->mutations[1]->pos,1.25);
  //BOOST_CHECK_EQUAL(offspring[2].second->mutations[1]->pos,1.5);
}

