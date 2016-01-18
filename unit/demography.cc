/*! \file demography.cc
  \ingroup unit 
  \brief Testing functions in fwdpp/demography.hpp and constructing a metapop from a singlepop
*/
#define BOOST_TEST_MODULE demography
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/metapop.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>

using mtype = KTfwd::popgenmut;
using singlepop_t = KTfwd::singlepop<mtype>;
using metapop_t = KTfwd::metapop<mtype>;

BOOST_AUTO_TEST_CASE( copy_construct_metapop_from_singlepop )
{
  singlepop_t pop(1000);
  BOOST_REQUIRE_EQUAL(pop.N,1000);
  metapop_t mpop(pop);
  BOOST_REQUIRE_EQUAL(mpop.Ns.size(),1);
  BOOST_REQUIRE_EQUAL(mpop.Ns[0],1000);
  BOOST_REQUIRE_EQUAL(mpop.diploids.size(),1);
  BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(),1000);
  BOOST_REQUIRE_EQUAL(pop.mutations==mpop.mutations,true);
  BOOST_REQUIRE_EQUAL(pop.gametes==mpop.gametes,true);
  BOOST_REQUIRE_EQUAL(pop.diploids==mpop.diploids[0],true);
}

BOOST_AUTO_TEST_CASE( move_construct_metapop_from_singlepop )
{
  singlepop_t pop(1000);
  BOOST_REQUIRE_EQUAL(pop.N,1000);
  metapop_t mpop(std::move(pop));
  BOOST_REQUIRE_EQUAL(mpop.Ns.size(),1);
  BOOST_REQUIRE_EQUAL(mpop.Ns[0],1000);
  BOOST_REQUIRE_EQUAL(mpop.diploids.size(),1);
  BOOST_REQUIRE_EQUAL(mpop.diploids[0].size(),1000);

  //We do not test sizes of elements in pop b/c
  //what happens will be a little bit compiler-dependent.
  //Typically, though, most or all of the elements in pop
  //should be empty.
}
