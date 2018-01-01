/*!
  \file sugar_singlepop.cc
  \ingroup unit
  \brief Testing fwdpp::singlepop
*/
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/io/serialize_population.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"

BOOST_FIXTURE_TEST_SUITE(test_singlepop, singlepop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(singlepop_sugar_test1)
{
    simulate_singlepop(pop);

    auto pop2(pop);

    BOOST_CHECK_EQUAL(pop == pop2, true);
}

BOOST_AUTO_TEST_CASE(singlepop_sugar_test3)
{
    simulate_singlepop(pop);

    auto pop2(std::move(pop));
    // Should be false b/c move will leave
    // pop's containers in a wacky state
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

BOOST_AUTO_TEST_CASE(singlepop_sugar_test4)
{
    simulate_singlepop(pop);

    auto pop2 = std::move(pop);
    // Should be false b/c move will leave
    // pop's containers in a wacky state
    BOOST_CHECK_EQUAL(pop == pop2, false);
}

//Test ability to serialize at different popsizes

BOOST_AUTO_TEST_CASE(singlepop_serialize_smallN)
{
	simulate_singlepop(pop,1000,1000);
    std::stringstream buffer;
	fwdpp::io::serialize_population(buffer, pop);
	singlepop_popgenmut_fixture::poptype pop2(0);
	fwdpp::io::deserialize_population(pop2,buffer);
	BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE(singlepop_serialize_mediumN)
{
	simulate_singlepop(pop,5000,5000);
    std::stringstream buffer;
	fwdpp::io::serialize_population(buffer, pop);
	singlepop_popgenmut_fixture::poptype pop2(0);
	fwdpp::io::deserialize_population(pop2,buffer);
	BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_CASE(singlepop_serialize_largeN)
{
	simulate_singlepop(pop,10000,10000);
    std::stringstream buffer;
	fwdpp::io::serialize_population(buffer, pop);
	singlepop_popgenmut_fixture::poptype pop2(0);
	fwdpp::io::deserialize_population(pop2,buffer);
	BOOST_CHECK_EQUAL(pop==pop2,true);
}

BOOST_AUTO_TEST_SUITE_END()
