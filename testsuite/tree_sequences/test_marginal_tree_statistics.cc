#include <iostream>
#include <fwdpp/ts/marginal_tree_functions/statistics.hpp>
#include <fwdpp/ts/decapitate.hpp>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection.hpp"
#include "simple_table_collection_polytomy.hpp"
#include "independent_implementations.hpp"

BOOST_AUTO_TEST_SUITE(test_marginal_tree_statistics)

BOOST_FIXTURE_TEST_CASE(test_total_time_simple, simple_table_collection)
{
    auto ttime = fwdpp::ts::total_time(tv.tree(), tables.nodes, true);
    BOOST_REQUIRE_EQUAL(ttime, total_time);
    BOOST_REQUIRE_EQUAL(
        ttime, naive_branch_length(tv.tree(), tables.nodes, true));
}

BOOST_FIXTURE_TEST_CASE(test_total_time_polytomy,
                        simple_table_collection_polytomy)
{
    auto ttime = fwdpp::ts::total_time(tv.tree(), tables.nodes, true);
    BOOST_REQUIRE_EQUAL(ttime, total_time);
}

BOOST_FIXTURE_TEST_CASE(test_total_time_polytomy_decapitate,
                        simple_table_collection_polytomy)
{
    fwdpp::ts::decapitate(tables, 0., false);
    reset_visitor(false);
    auto ttime = fwdpp::ts::total_time(tv.tree(), tables.nodes, true);
    BOOST_REQUIRE(ttime != total_time);
    BOOST_REQUIRE_EQUAL(
        ttime, naive_branch_length(tv.tree(), tables.nodes, true));
}

BOOST_AUTO_TEST_SUITE_END()

