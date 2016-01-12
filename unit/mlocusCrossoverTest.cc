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
#include <vector>
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
using gvector = std::vector<gtype>;
using mutvector = std::vector<mut>;
using diploid_t = std::vector< std::pair< std::size_t,std::size_t> >;

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
void setup3locus2( gvector & gametes,
		   mutvector & mvector,
		   diploid_t & diploid )
{
  gametes = gvector(6,gtype(1));
  mvector.clear();

  //To set this up, let's add the mutations:
  mvector.emplace_back(0.5,0.);
  mvector.emplace_back(0.75,0.);
  mvector.emplace_back(0.25,0.);
  mvector.emplace_back(0.9,0.);
  mvector.emplace_back(1.25,0.1);
  mvector.emplace_back(1.5,0.1);
  mvector.emplace_back(1.1,0.1);
 
  //put the mutations into gametes
  gametes[0].mutations.push_back(0);  //gamete 1, locus 1
  gametes[1].mutations.push_back(2);    //gamete 2, locus 1
  gametes[2].mutations.push_back(1);
  gametes[3].mutations.push_back(3);
  gametes[4].mutations.push_back(4);
  gametes[4].mutations.push_back(5);
  gametes[5].mutations.push_back(6);

  //Now, make a diploid
  diploid.clear();
  diploid.emplace_back(std::make_pair(0,1));
  diploid.emplace_back(std::make_pair(2,3));
  diploid.emplace_back(std::make_pair(4,5));
  return;
}

BOOST_AUTO_TEST_CASE( three_locus_test_1 )
{
  gvector gametes;
  mutvector mvector;
  diploid_t diploid;            //parent 1
  setup3locus2(gametes,mvector,diploid);
  std::vector<unsigned> mcounts(mvector.size(),1);
  diploid_t diploid2(diploid); //parent 2

  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  std::vector<std::vector<double> > rec1 { std::vector<double>{ 0.3,MVAL },
      std::vector<double>{ 0.55,0.8,MVAL },
	std::vector<double>{1.3,MVAL}
  };

  //We use these to "fake" what we want to happen between loci.
  std::vector<double> r_bw_loci = {1.,0.};
  std::vector<diploid_t> diploids({diploid});

  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(gametes,mvector);
  auto mutation_recycling_bin = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
  auto gamete_recycling_bin = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3

  std::vector<std::function<std::vector<double>(const gtype &,
						const gtype &,
						const mutvector &)> > recpols
  {
    [&rec1](const gtype &,
	    const gtype &,
	    const mutvector &) { return rec1[0]; },
      [&rec1](const gtype &,
	      const gtype &,
	      const mutvector &) { return rec1[1]; },
	[&rec1](const gtype &,
		const gtype &,
		const mutvector &) { return rec1[2]; }
  };

  auto offspring  = KTfwd::fwdpp_internal::multilocus_rec_mut(r,
							      diploid,diploid2,
							      mutation_recycling_bin,gamete_recycling_bin,gamete_lookup,
							      recpols,
							      [](gsl_rng * __r, const double & __d) { return __d; },
							      &r_bw_loci[0],
							      0,0,gametes,mvector,neutral,selected);
 
  BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(),0);
  BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(),0);
  BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(),1);
  BOOST_CHECK_EQUAL(mvector[gametes[offspring[2].first].mutations[0]].pos,1.25);
}

BOOST_AUTO_TEST_CASE( three_locus_test_2 )
{
  gvector gametes;
  mutvector mvector;
  diploid_t diploid;            //parent 1
  setup3locus2(gametes,mvector,diploid);
  std::vector<unsigned> mcounts(mvector.size(),1);
  diploid_t diploid2(diploid); //parent 2

  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
   auto rec1 = std::vector< std::vector<double> > { std::vector<double>{ 0.45, MVAL },
						    std::vector<double>{ MVAL },
						    std::vector<double>{ MVAL }
   };

  //We use these to "fake" what we want to happen between loci.
  std::vector<double> r_bw_loci = {1.,0.};
  std::vector<diploid_t> diploids({diploid});

  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(gametes,mvector);
  auto mutation_recycling_bin = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
  auto gamete_recycling_bin = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3

  std::vector<std::function<std::vector<double>(const gtype &,
						const gtype &,
						const mutvector &)> > recpols
  {
    [&rec1](const gtype &,
	    const gtype &,
	    const mutvector &) { return rec1[0]; },
      [&rec1](const gtype &,
	      const gtype &,
	      const mutvector &) { return rec1[1]; },
	[&rec1](const gtype &,
		const gtype &,
		const mutvector &) { return rec1[2]; }
  };

  auto offspring  = KTfwd::fwdpp_internal::multilocus_rec_mut(r,
							      diploid,diploid2,
							      mutation_recycling_bin,gamete_recycling_bin,gamete_lookup,
							      recpols,
							      [](gsl_rng * __r, const double & __d) { return __d; },
							      &r_bw_loci[0],
							      0,0,gametes,mvector,neutral,selected);
  
  BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(),0);
  BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(),1);
  BOOST_CHECK_EQUAL(mvector[gametes[offspring[1].first].mutations[0]].pos,0.75);
  BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(),2);
  BOOST_CHECK_EQUAL(mvector[gametes[offspring[2].first].mutations[0]].pos,1.25);
  BOOST_CHECK_EQUAL(mvector[gametes[offspring[2].first].mutations[1]].pos,1.5);
 
}

BOOST_AUTO_TEST_CASE( three_locus_test_3 )
{
  gvector gametes;
  mutvector mvector;
  diploid_t diploid;            //parent 1
  setup3locus2(gametes,mvector,diploid);
  std::vector<unsigned> mcounts(mvector.size(),1);
  diploid_t diploid2(diploid); //parent 2

  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector< std::vector<double> > { std::vector<double>{ 0.1,0.45, MVAL },
 						   std::vector<double>{ MVAL },
 						   std::vector<double>{ MVAL }
  };

  //We use these to "fake" what we want to happen between loci.
  std::vector<double> r_bw_loci = {1.,1.};
  std::vector<diploid_t> diploids({diploid});

  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(gametes,mvector);
  auto mutation_recycling_bin = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
  auto gamete_recycling_bin = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3

  std::vector<std::function<std::vector<double>(const gtype &,
						const gtype &,
						const mutvector &)> > recpols
  {
    [&rec1](const gtype &,
	    const gtype &,
	    const mutvector &) { return rec1[0]; },
      [&rec1](const gtype &,
	      const gtype &,
	      const mutvector &) { return rec1[1]; },
	[&rec1](const gtype &,
		const gtype &,
		const mutvector &) { return rec1[2]; }
  };

  auto offspring  = KTfwd::fwdpp_internal::multilocus_rec_mut(r,
							      diploid,diploid2,
							      mutation_recycling_bin,gamete_recycling_bin,gamete_lookup,
							      recpols,
							      [](gsl_rng * __r, const double & __d) { return __d; },
							      &r_bw_loci[0],
							      0,0,gametes,mvector,neutral,selected);

  BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(),2);
  BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(),1);
  BOOST_CHECK_EQUAL(mvector[gametes[offspring[1].first].mutations[0]].pos,0.9);
  BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(),2);
  BOOST_CHECK_EQUAL(mvector[gametes[offspring[2].first].mutations[0]].pos,1.25);
  BOOST_CHECK_EQUAL(mvector[gametes[offspring[2].first].mutations[1]].pos,1.5);
  
}

BOOST_AUTO_TEST_CASE( three_locus_test_4 )
{
  gvector gametes;
  mutvector mvector;
  diploid_t diploid;            //parent 1
  setup3locus2(gametes,mvector,diploid);
  std::vector<unsigned> mcounts(mvector.size(),1);
  diploid_t diploid2(diploid); //parent 2

  //positions of x-overs within loci
  auto MVAL=std::numeric_limits<double>::max();
  auto rec1 = std::vector< std::vector<double> > { std::vector<double>{ 0.1,0.45, 0.51, MVAL },
 						   std::vector<double>{ 0.6,0.7,0.99, MVAL },
 						   std::vector<double>{ 1,1.2,1.3,1.55,MVAL }
  };

  //We use these to "fake" what we want to happen between loci.
  std::vector<double> r_bw_loci = {1.,1.};
  std::vector<diploid_t> diploids({diploid});

  auto gamete_lookup = KTfwd::fwdpp_internal::gamete_lookup_table(gametes,mvector);
  auto mutation_recycling_bin = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
  auto gamete_recycling_bin = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
  gtype::mutation_container neutral,selected; //req'd as of 0.3.3

  std::vector<std::function<std::vector<double>(const gtype &,
						const gtype &,
						const mutvector &)> > recpols
  {
    [&rec1](const gtype &,
	    const gtype &,
	    const mutvector &) { return rec1[0]; },
      [&rec1](const gtype &,
	      const gtype &,
	      const mutvector &) { return rec1[1]; },
	[&rec1](const gtype &,
		const gtype &,
		const mutvector &) { return rec1[2]; }
  };

  auto offspring  = KTfwd::fwdpp_internal::multilocus_rec_mut(r,
							      diploid,diploid2,
							      mutation_recycling_bin,gamete_recycling_bin,gamete_lookup,
							      recpols,
							      [](gsl_rng * __r, const double & __d) { return __d; },
							      &r_bw_loci[0],
							      0,0,gametes,mvector,neutral,selected);
       BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(),2);
       BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(),1);
       BOOST_CHECK_EQUAL(mvector[gametes[offspring[1].first].mutations[0]].pos,0.75);
       BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(),2);
       BOOST_CHECK_EQUAL(mvector[gametes[offspring[2].first].mutations[0]].pos,1.1);
       BOOST_CHECK_EQUAL(mvector[gametes[offspring[2].first].mutations[1]].pos,1.25);
}
