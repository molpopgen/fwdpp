#include <stdexcept>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/marginal_tree_functions/children.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include "simple_table_collection.hpp"
#include "independent_implementations.hpp"

BOOST_AUTO_TEST_SUITE(test_sample_traversal)

BOOST_FIXTURE_TEST_CASE(num_samples_simple_table_collection,
                        simple_table_collection)
{
    // NOTE: only works b/c no recombination, so node table size
    // is the number of nodes in the tree.  More generally,
    // we need to do a node traveral, but that hasn't been tested yet
    // (see below).
    for (std::size_t i = 0; i < tables.node_table.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(fwdpp::ts::num_samples(tv.tree(), i),
                                naive_num_samples(tv.tree(), i));
        }
}

BOOST_FIXTURE_TEST_CASE(test_sample_lists, simple_table_collection)
{
    for (std::size_t i = 0; i < tables.node_table.size(); ++i)
        {
            auto x = fwdpp::ts::get_samples(tv.tree(), i);
            auto y = naive_get_samples(tv.tree(), i);
            BOOST_CHECK_EQUAL(x.size(), y.size());
            BOOST_CHECK(fwdpp::ts::get_samples(tv.tree(), i)
                        == naive_get_samples(tv.tree(), i));
        }
}

BOOST_FIXTURE_TEST_CASE(test_process_samples_as_nodes, simple_table_collection)
{
    samples = { 1, 2 };
    reset_visitor(true);
    std::vector<fwdpp::ts::TS_NODE_INT> s;
    fwdpp::ts::process_samples(
        tv.tree(), fwdpp::ts::convert_sample_index_to_nodes(true), 6,
        [&s](fwdpp::ts::TS_NODE_INT u) { s.push_back(u); });
    BOOST_REQUIRE(s == samples);
}

BOOST_FIXTURE_TEST_CASE(test_process_samples_as_indexes,
                        simple_table_collection)
{
    samples = { 1, 2 };
    reset_visitor(true);
    std::vector<fwdpp::ts::TS_NODE_INT> s;
    fwdpp::ts::process_samples(
        tv.tree(), fwdpp::ts::convert_sample_index_to_nodes(false), 6,
        [&s](fwdpp::ts::TS_NODE_INT u) { s.push_back(u); });
    for (auto& i : samples)
        {
            i -= 1;
        }
    BOOST_REQUIRE(s == samples);
}

BOOST_FIXTURE_TEST_CASE(test_subset_of_nodes_as_samples,
                        simple_table_collection)
{
    samples = { 0, 1, 3 };
    reset_visitor_and_samples(samples, true);
    BOOST_REQUIRE_EQUAL(naive_num_samples(tv.tree(), 6), samples.size());
    auto x = fwdpp::ts::get_samples(tv.tree(), 6);
    BOOST_CHECK(x.size() == samples.size());
    BOOST_CHECK(x == samples);
}

BOOST_FIXTURE_TEST_CASE(test_exception, simple_table_collection)
{
    // need a visitor that does NOT update samples lists
    BOOST_REQUIRE_THROW(
        {
            reset_visitor(false);
            // Throws b/c tv is not updating samples lists
            fwdpp::ts::num_samples(tv.tree(), 6);
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
