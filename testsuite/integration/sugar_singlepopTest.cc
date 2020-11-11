/*!
  \file sugar_singlepopTest.cc
  \ingroup unit
  \brief Testing fwdpp::diploid_population
*/
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/forward_types_serialization.hpp>
#include <fwdpp/io/serialize_population.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"

BOOST_FIXTURE_TEST_SUITE(test_diploid_population, diploid_population_mutation_fixture)

BOOST_AUTO_TEST_CASE(diploid_population_sugar_test1)
{
    simulate_diploid_population(pop);

    auto pop2(pop);

    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(diploid_population_sugar_test3)
{
    simulate_diploid_population(pop);

    auto pop2(std::move(pop));
    // Should be false b/c move will leave
    // pop's containers in a wacky state
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(diploid_population_sugar_test4)
{
    simulate_diploid_population(pop);

    auto pop2 = std::move(pop);
    // Should be false b/c move will leave
    // pop's containers in a wacky state
    BOOST_CHECK_EQUAL(pop == pop2, false);
}


BOOST_AUTO_TEST_CASE(diploid_population_serialize)
{
    simulate_diploid_population(pop, 1000, 1000);
    std::stringstream buffer;
    fwdpp::io::serialize_population(buffer, pop);
    diploid_population_mutation_fixture::poptype pop2(0);
    fwdpp::io::deserialize_population(buffer, pop2);
    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_SUITE_END()
