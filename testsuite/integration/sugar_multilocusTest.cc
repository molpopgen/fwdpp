/*!
  \file sugar_multilocus.cc
  \ingroup unit
  \brief Testing fwdpp::multiloc
*/
#include <unistd.h>
#include <config.h>
#include <zlib.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/multiloc.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"

using poptype = multiloc_popgenmut_fixture::poptype;

BOOST_AUTO_TEST_CASE(multiloc_sugar_test1)
{
    multiloc_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       multiloc_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2(f.pop);
    BOOST_CHECK_EQUAL(f.pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(multiloc_sugar_test2)
{
    multiloc_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       multiloc_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2(0, 0);
    std::stringstream buffer;
	fwdpp::serialize_population(buffer,f.pop);
    fwdpp::deserialize_population(pop2, buffer);
    BOOST_CHECK_EQUAL(f.pop == pop2, true);
}



BOOST_AUTO_TEST_CASE(multiloc_sugar_test3)
{
    multiloc_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       multiloc_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2(std::move(f.pop));
    BOOST_CHECK_EQUAL(f.pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(multiloc_sugar_test4)
{
    multiloc_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       multiloc_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2 = std::move(f.pop);
    BOOST_CHECK_EQUAL(f.pop == pop2, false);
}
