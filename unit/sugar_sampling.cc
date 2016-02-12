/*! 
  \file sugar_sampling.cc 
  \ingroup unit 
  \brief Testing KTfwd::sample and KTfwd::sample_separate
*/
#define BOOST_TEST_MODULE sugar_sampling
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mutation_t = KTfwd::popgenmut;

using singlepop_t = KTfwd::singlepop<mutation_t>;
using metapop_t = KTfwd::metapop<mutation_t>;
using multiloc_t = KTfwd::multiloc<mutation_t>;
  
KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

//This is an API check
BOOST_AUTO_TEST_CASE( singledeme1 )
{
  singlepop_t pop(1000);
  auto s = KTfwd::sample_separate(rng.get(),pop,10,true);
}

//This is an API check
BOOST_AUTO_TEST_CASE( singledeme2 )
{
  singlepop_t pop(1000);
  auto s = KTfwd::sample_separate(pop,std::vector<unsigned>({0,1,2,999}),true);
}

BOOST_AUTO_TEST_CASE( singledeme3 )
{
  singlepop_t pop(1000);
  auto s = KTfwd::sample_separate(pop,std::vector<unsigned>(),true);
  BOOST_REQUIRE(s.first.empty()==true);
  BOOST_REQUIRE(s.second.empty()==true);
}

BOOST_AUTO_TEST_CASE( singledeme4 )
{
  singlepop_t pop(1000);
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(pop,std::vector<unsigned>({1000}),true),
		      std::out_of_range);
}

//This is an API check
BOOST_AUTO_TEST_CASE( singledeme5 )
{
  singlepop_t pop(1000);
  auto s = KTfwd::sample(rng.get(),pop,10,true);
}

//This is an API check
BOOST_AUTO_TEST_CASE( singledeme6 )
{
  singlepop_t pop(1000);
  auto s = KTfwd::sample(pop,std::vector<unsigned>({0,1,2,999}),true);
}

BOOST_AUTO_TEST_CASE( singledeme7 )
{
  singlepop_t pop(1000);
  auto s = KTfwd::sample(pop,std::vector<unsigned>(),true);
  BOOST_REQUIRE(s.empty()==true);
}

BOOST_AUTO_TEST_CASE( singledeme8 )
{
  singlepop_t pop(1000);
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample(pop,std::vector<unsigned>({1000}),true),
		      std::out_of_range);
}

//This is an API check
BOOST_AUTO_TEST_CASE( metapop1 )
{
  metapop_t pop({1000,1000});
  auto s = KTfwd::sample_separate(rng.get(),pop,0,10,true);
}
BOOST_AUTO_TEST_CASE( metapop1_throw )
{
  metapop_t pop({1000,1000});
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(rng.get(),pop,2,10,true),
		      std::out_of_range);
}

//This is an API check
BOOST_AUTO_TEST_CASE( metapop2 )
{
  metapop_t pop({1000,1000});
  auto s = KTfwd::sample_separate(pop,0,std::vector<unsigned>({10,20,50}),true);
}

BOOST_AUTO_TEST_CASE( metapop2_throw )
{
  metapop_t pop({1000,1000});
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(pop,2,std::vector<unsigned>({10,20,50}),true),
		      std::out_of_range);
}

BOOST_AUTO_TEST_CASE( metapop2_throw2 )
{
  metapop_t pop({1000,1000});
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(pop,0,std::vector<unsigned>({10,20,50,1000}),true),
		      std::out_of_range);
}

 BOOST_AUTO_TEST_CASE( metapop3 )
 {
   metapop_t pop({1000,1000});
   auto s = KTfwd::sample_separate(pop,0,std::vector<unsigned>(),true);
   BOOST_REQUIRE(s.first.empty()==true);
   BOOST_REQUIRE(s.second.empty()==true);
}

BOOST_AUTO_TEST_CASE( metapop4 )
{
  metapop_t pop({1000,1000});
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample_separate(pop,2,std::vector<unsigned>({1000}),true),
		      std::out_of_range);
}

//This is an API check
BOOST_AUTO_TEST_CASE( metapop5 )
{
  metapop_t pop({1000,1000});
  auto s = KTfwd::sample(rng.get(),pop,0,10,true);
}

BOOST_AUTO_TEST_CASE( metapop5_throw )
{
  metapop_t pop({1000,1000});
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample(rng.get(),pop,2,10,true),
		      std::out_of_range);
}

//This is an API check
BOOST_AUTO_TEST_CASE( metapop6 )
{
  metapop_t pop({1000,1000});
  auto s = KTfwd::sample(pop,0,std::vector<unsigned>({0,1,2,999}),true);
}

BOOST_AUTO_TEST_CASE( metapop7 )
{
  metapop_t pop({1000,1000});
  auto s = KTfwd::sample(pop,0,std::vector<unsigned>(),true);
  BOOST_REQUIRE(s.empty()==true);
}

BOOST_AUTO_TEST_CASE( metapop8 )
{
  metapop_t pop({1000,1000});
  BOOST_REQUIRE_THROW(auto s = KTfwd::sample(pop,2,std::vector<unsigned>({1000}),true),
		      std::out_of_range);
}

//Correctness tests

BOOST_AUTO_TEST_CASE( singlepop_1 )
{
  singlepop_t pop(1000);
  pop.mutations.emplace_back( 0.1,0,0,0 );
  pop.mcounts.emplace_back(1);
  pop.gametes.emplace_back( 1, std::vector<std::size_t>{0}, std::vector<std::size_t>{} );
  pop.gametes[0].n--;
  BOOST_REQUIRE_EQUAL( KTfwd::check_sum(pop.gametes,2000), true );
  BOOST_REQUIRE_EQUAL( KTfwd::popdata_sane(pop.diploids,pop.gametes,pop.mcounts), true );
  pop.diploids[1].first=1;

  auto x = KTfwd::sample(pop,{0,1},true);

  BOOST_REQUIRE_EQUAL(x.size(),1);
  BOOST_REQUIRE_EQUAL(std::count(x[0].second.begin(),x[0].second.end(),'1'),1);
  BOOST_REQUIRE_EQUAL(x[0].second[2],'1');
  //Now, make that diploid homozygous for the mutation
  /*
    FROM HERE ON OUT, INTERNAL DATA STRUCTURES ARE INCONSISTENT.
    This doesn't matter for these tests, but it could in future releases, e.g.,
    if internal code adds some assertions...
  */
  pop.diploids[1].second=1;

  x = KTfwd::sample(pop,{0,1},true);
  BOOST_REQUIRE_EQUAL(x.size(),1);
  BOOST_REQUIRE_EQUAL(std::count(x[0].second.begin(),x[0].second.end(),'1'),2);
  BOOST_REQUIRE_EQUAL(x[0].second[2],'1');
  BOOST_REQUIRE_EQUAL(x[0].second[3],'1');

  //Now, we are only sampling individual 1, and thus
  //the variant at 0.1 will be fixed in the sample,
  //and the sample should be empty
  x = KTfwd::sample(pop,{1},true);
  BOOST_REQUIRE_EQUAL(x.empty(),true);

  //Now, allow fixed variants in the sample
  x = KTfwd::sample(pop,{1},false);
  BOOST_REQUIRE_EQUAL(x.empty(),false);

  //now, add a fixation with position < 0.1
  pop.fixations.emplace_back(-0.1,0,0,0);

  //Now, we sample, but we don't want the fixations...
  x = KTfwd::sample(pop,{0,1},true);
  BOOST_REQUIRE_EQUAL(x.size(),1);
  BOOST_REQUIRE_EQUAL(std::count(x[0].second.begin(),x[0].second.end(),'1'),2);
  BOOST_REQUIRE_EQUAL(x[0].second[2],'1');
  BOOST_REQUIRE_EQUAL(x[0].second[3],'1');

  //Now, we want the fixation
  x = KTfwd::sample(pop,{0,1},false);
  BOOST_REQUIRE_EQUAL(x.size(),2);
  BOOST_REQUIRE_EQUAL(x[0].first,-0.1);
  BOOST_REQUIRE_EQUAL(x[1].first,0.1);
  BOOST_REQUIRE_EQUAL( std::is_sorted(x.begin(),x.end(),[](const KTfwd::sample_site_t & i,const KTfwd::sample_site_t & j) {
	return i.first<j.first;
      }), true);
}


//Now, a metapop
BOOST_AUTO_TEST_CASE( metapop_1 )
{
  metapop_t mpop({1000,1000});
  mpop.fixations.emplace_back( -0.1,0,0,0 );

  //Make sure that fiations end up in samples
  auto x = KTfwd::sample_separate( mpop, 0, {0,1,2}, false );
  BOOST_REQUIRE_EQUAL(x.first.size(),1);

  //Add a selected fixation
  mpop.fixations.emplace_back( -0.2,1.0,0,0 );
  x = KTfwd::sample_separate( mpop, 0, {0,1,2}, false );
  BOOST_REQUIRE_EQUAL(x.first.size(),1);
}
