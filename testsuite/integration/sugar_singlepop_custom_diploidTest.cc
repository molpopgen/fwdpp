/*!
  \file sugar_singlepop_custom_diploidTest.cc
  \ingroup unit
  \brief Testing single-deme sugar functionality with custom diploids
*/
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>
#include <testsuite/util/quick_evolve_sugar.hpp>

using mutation_t = fwdpp::popgenmut;

BOOST_FIXTURE_TEST_SUITE(test_diploid_population_custom,
                         diploid_population_popgenmut_custom_fixture)

BOOST_AUTO_TEST_CASE(diploid_population_sugar_custom_test1)
{
    simulate_diploid_population(pop);

    auto pop2(pop);

    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(diploid_population_sugar_custom_test3)
{
    simulate_diploid_population(pop);

    auto pop2(std::move(pop));
    // Should be false b/c move will leave
    // pop's containers in a wacky state
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(diploid_population_sugar_custom_test4)
{
    simulate_diploid_population(pop);

    auto pop2 = std::move(pop);
    // Should be false b/c move will leave
    // pop's containers in a wacky state
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_SUITE_END()
