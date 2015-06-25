/*!
  \file mlocusCrossoverTest.cc
  \ingroup unit
  \brief Tests KTfwd::fwdpp_internal::multilocus_rec
*/

#define BOOST_TEST_MODULE mlocusCrossoverTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <iostream>
//For this unit test, this symbol eliminates the mutation-related part of KTfwd::fwdpp_internal::multiloc_rec_mut,
//which means we don't have to write as much boilerplate code to test the more complex logic.
//Plus, mutation stuff is unit-tested elsewhere
#define FWDPP_UNIT_TESTING
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
void setup3locus2( glist & gametes,
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
  setup3locus2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2
  // //This block makes sure that setup3locus2 is working as far as gametes/mutations:
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

  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  std::vector<std::vector<double> > rec1 { std::vector<double>{ 0.3,MVAL },
      std::vector<double>{ 0.55,0.8,MVAL },
	std::vector<double>{1.3,MVAL}
  };

  //We use these to "fake" what we want to happen between loci.
  std::vector<double> r_bw_loci = {1.,0.};
  std::vector<diploid_t> diploids({diploid});
  auto offspring = diploids.begin(); 
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3

  std::function<unsigned(glist::iterator &, glist::iterator &, decltype(gamete_lookup) &)> recpol1 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
    return KTfwd::recombine_gametes(rec1[0],&gametes,g1,g2,lookup,neutral,selected);
  },
    recpol2 =  [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
      return KTfwd::recombine_gametes(rec1[1],&gametes,g1,g2,lookup,neutral,selected);
    },
      recpol3 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
	return KTfwd::recombine_gametes(rec1[2],&gametes,g1,g2,lookup,neutral,selected);
      };
      
      std::vector< decltype(recpol1) >  recpols ({ recpol1,recpol2,recpol3 });
  KTfwd::fwdpp_internal::multilocus_rec_mut(r,
					    diploid,diploid2,offspring,gamete_lookup,
					    recpols,
					    [](gsl_rng * __r, const double & __d) { return __d; },
					    &r_bw_loci[0],
					    0,0);
 
  BOOST_CHECK_EQUAL((*offspring)[0].first->mutations.size(),0);
  BOOST_CHECK_EQUAL((*offspring)[1].first->mutations.size(),0);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations.size(),1);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations[0]->pos,1.25);
}

BOOST_AUTO_TEST_CASE( three_locus_test_2 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;            //parent 1
  setup3locus2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2


  /*
    Make vectors of events regarding what happens withn and between each
    locus, for each parent.
  */
  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();

  auto rec1 = std::vector< std::vector<double> > { std::vector<double>{ 0.45, MVAL },
						   std::vector<double>{ MVAL },
						   std::vector<double>{ MVAL }
  };
 
 //We use these to "fake" what we want to happen between loci.
  std::vector<double> r_bw_loci = {1.,0.};
  std::vector<diploid_t> diploids({diploid});
  auto offspring = diploids.begin(); 
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  
  std::function<unsigned(glist::iterator &, glist::iterator &, decltype(gamete_lookup) &)> recpol1 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
    return KTfwd::recombine_gametes(rec1[0],&gametes,g1,g2,lookup,neutral,selected);
  },
    recpol2 =  [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
      return KTfwd::recombine_gametes(rec1[1],&gametes,g1,g2,lookup,neutral,selected);
    },
      recpol3 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
	return KTfwd::recombine_gametes(rec1[2],&gametes,g1,g2,lookup,neutral,selected);
      };
      
      std::vector< decltype(recpol1) >  recpols ({ recpol1,recpol2,recpol3 });
      KTfwd::fwdpp_internal::multilocus_rec_mut(r,
						diploid,diploid2,offspring,gamete_lookup,
						recpols,
						[](gsl_rng * __r, const double & __d) { return __d; },
						&r_bw_loci[0],
						0,0);

  BOOST_CHECK_EQUAL((*offspring)[0].first->mutations.size(),0);
  BOOST_CHECK_EQUAL((*offspring)[1].first->mutations.size(),1);
  BOOST_CHECK_EQUAL((*offspring)[1].first->mutations[0]->pos,0.75);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations.size(),2);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations[0]->pos,1.25);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations[1]->pos,1.5);
  
 
}

BOOST_AUTO_TEST_CASE( three_locus_test_3 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;            //parent 1
  setup3locus2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2


  /*
    Make vectors of events regarding what happens withn and between each
    locus, for each parent.
  */
  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector< std::vector<double> > { std::vector<double>{ 0.1,0.45, MVAL },
						   std::vector<double>{ MVAL },
						   std::vector<double>{ MVAL }
  };
  //recombinations b/w all loci this time...
  std::vector<double> r_bw_loci = {1.,1.};
  std::vector<diploid_t> diploids({diploid});
  auto offspring = diploids.begin(); 
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  
  std::function<unsigned(glist::iterator &, glist::iterator &, decltype(gamete_lookup) &)> recpol1 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
    return KTfwd::recombine_gametes(rec1[0],&gametes,g1,g2,lookup,neutral,selected);
  },
    recpol2 =  [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
      return KTfwd::recombine_gametes(rec1[1],&gametes,g1,g2,lookup,neutral,selected);
    },
      recpol3 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
	return KTfwd::recombine_gametes(rec1[2],&gametes,g1,g2,lookup,neutral,selected);
      };
      
      std::vector< decltype(recpol1) >  recpols ({ recpol1,recpol2,recpol3 });
      KTfwd::fwdpp_internal::multilocus_rec_mut(r,
						diploid,diploid2,offspring,gamete_lookup,
						recpols,
						[](gsl_rng * __r, const double & __d) { return __d; },
						&r_bw_loci[0],
						0,0);

  
  BOOST_CHECK_EQUAL((*offspring)[0].first->mutations.size(),2);
  BOOST_CHECK_EQUAL((*offspring)[1].first->mutations.size(),1);
  BOOST_CHECK_EQUAL((*offspring)[1].first->mutations[0]->pos,0.9);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations.size(),2);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations[0]->pos,1.25);
  BOOST_CHECK_EQUAL((*offspring)[2].first->mutations[1]->pos,1.5);
}

BOOST_AUTO_TEST_CASE( three_locus_test_4 )
{
  glist gametes;
  mutlist mlist;
  diploid_t diploid;            //parent 1
  setup3locus2(gametes,mlist,diploid);
  diploid_t diploid2(diploid); //parent 2

  /*
    Make vectors of events regarding what happens withn and between each
    locus, for each parent.
  */
  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector< std::vector<double> > { std::vector<double>{ 0.1,0.45, 0.51, MVAL },
						   std::vector<double>{ 0.6,0.7,0.99, MVAL },
						   std::vector<double>{ 1,1.2,1.3,1.55,MVAL }
  };
  //recombinations b/w all loci this time...
  std::vector<double> r_bw_loci = {1.,1.};
   std::vector<diploid_t> diploids({diploid});
  auto offspring = diploids.begin(); 
  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(&gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3
  
  std::function<unsigned(glist::iterator &, glist::iterator &, decltype(gamete_lookup) &)> recpol1 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
    return KTfwd::recombine_gametes(rec1[0],&gametes,g1,g2,lookup,neutral,selected);
  },
    recpol2 =  [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
      return KTfwd::recombine_gametes(rec1[1],&gametes,g1,g2,lookup,neutral,selected);
    },
      recpol3 = [&rec1,&gametes,&neutral,&selected](glist::iterator & g1, glist::iterator & g2, decltype(gamete_lookup) & lookup) {
	return KTfwd::recombine_gametes(rec1[2],&gametes,g1,g2,lookup,neutral,selected);
      };
      
      std::vector< decltype(recpol1) >  recpols ({ recpol1,recpol2,recpol3 });
      KTfwd::fwdpp_internal::multilocus_rec_mut(r,
						diploid,diploid2,offspring,gamete_lookup,
						recpols,
						[](gsl_rng * __r, const double & __d) { return __d; },
						&r_bw_loci[0],
						0,0);
  
  BOOST_CHECK_EQUAL((*(offspring))[0].first->mutations.size(),2);
  BOOST_CHECK_EQUAL((*(offspring))[1].first->mutations.size(),1);
  BOOST_CHECK_EQUAL((*(offspring))[1].first->mutations[0]->pos,0.75);
  BOOST_CHECK_EQUAL((*(offspring))[2].first->mutations.size(),2);
  BOOST_CHECK_EQUAL((*(offspring))[2].first->mutations[0]->pos,1.1);
  BOOST_CHECK_EQUAL((*(offspring))[2].first->mutations[1]->pos,1.25);
}

