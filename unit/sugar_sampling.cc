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
#include <fwdpp/sugar/sampling.hpp>
#include <fwdpp/sugar/infsites.hpp>

using mutation_t = KTfwd::popgenmut;
using mwriter = KTfwd::mutation_writer;
using mreader = KTfwd::mutation_reader<mutation_t>;

using singlepop_t = KTfwd::singlepop<mutation_t>;

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
