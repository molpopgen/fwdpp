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
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/mlocuspop.hpp>
#include <fwdpp/forward_types_serialization.hpp>
#include <fwdpp/io/serialize_population.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"

BOOST_FIXTURE_TEST_SUITE(mlocus_tests, mlocuspop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test1)
{
    simulate_mlocuspop(pop, rng, mutmodels, recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), mu,
                       rbw, generation);
    decltype(pop) pop2(pop);
    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test2)
{
    simulate_mlocuspop(pop, rng, mutmodels, recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), mu,
                       rbw, generation);
    decltype(pop) pop2(0, decltype(poptype::locus_boundaries)());
    std::stringstream buffer;
    fwdpp::io::serialize_population(buffer, pop);
    fwdpp::io::deserialize_population(buffer, pop2);
    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test3)
{
    simulate_mlocuspop(pop, rng, mutmodels, recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), mu,
                       rbw, generation);
    auto pop2(std::move(pop));
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(mlocuspop_sugar_test4)
{
    simulate_mlocuspop(pop, rng, mutmodels, recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), mu,
                       rbw, generation);
    auto pop2 = std::move(pop);
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_SUITE_END()
