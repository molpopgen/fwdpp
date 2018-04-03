/*!
  \file sugar_multilocusTest.cc
  \ingroup unit
  \brief Testing fwdpp::mlocuspop
*/
#include <unistd.h>
#include <config.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/sugar/mlocuspop.hpp>
#include <fwdpp/forward_types_serialization.hpp>
#include <fwdpp/io/serialize_population.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"

using poptype = mlocuspop_popgenmut_fixture::poptype;

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test1)
{
    mlocuspop_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2(f.pop);
    BOOST_CHECK_EQUAL(f.pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test2)
{
    mlocuspop_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2(0, 0);
    std::stringstream buffer;
    fwdpp::io::serialize_population(buffer, f.pop);
    fwdpp::io::deserialize_population(buffer, pop2);
    BOOST_CHECK_EQUAL(f.pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test3)
{
    mlocuspop_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2(std::move(f.pop));
    BOOST_CHECK_EQUAL(f.pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test4)
{
    mlocuspop_popgenmut_fixture f;
    simulate_mlocuspop(f.pop, f.rng, f.mutmodels, f.recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), f.mu,
                       f.rbw, f.generation);
    poptype pop2 = std::move(f.pop);
    BOOST_CHECK_EQUAL(f.pop == pop2, false);
}
