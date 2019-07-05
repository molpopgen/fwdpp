#include <fwdpp/ts/tree_visitor.hpp>
#include <boost/test/unit_test.hpp>
#include "empty_table_collection.hpp"
#include "simple_table_collection.hpp"

BOOST_FIXTURE_TEST_SUITE(test_tree_visitor_with_simple_table_collection,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(test_construction)
{
    // NOTE: this is a "cheat", as we cannot get this far
    // if this doesn't work, as the fixture constructor
    // does this.
    BOOST_REQUIRE_NO_THROW({ fwdpp::ts::tree_visitor(tables, samples); });
}

BOOST_AUTO_TEST_CASE(test_unindexed_tables)
{
    tables.input_left.clear();
    tables.output_right.clear();
    BOOST_REQUIRE_THROW(fwdpp::ts::tree_visitor tv(tables, samples),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_tree_visitor_with_empty_table_collection,
                         empty_table_collection)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_THROW(fwdpp::ts::tree_visitor tv(tables, empty_samples),
                        fwdpp::ts::empty_samples);
}

BOOST_AUTO_TEST_SUITE_END()
