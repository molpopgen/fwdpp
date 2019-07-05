#include <iostream>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection.hpp"
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/marginal_tree_functions.hpp>

BOOST_FIXTURE_TEST_SUITE(test_simple_table_collection_traversal,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_NO_THROW({ fwdpp::ts::node_iterator(tv.tree()); });
}

BOOST_AUTO_TEST_CASE(test_preorder_traversal)
{
    fwdpp::ts::node_iterator ni(tv.tree());
    std::vector<fwdpp::ts::TS_NODE_INT> nodes;
    auto n = ni(fwdpp::ts::nodes_preorder());
    while (n != fwdpp::ts::TS_NULL_NODE)
        {
            nodes.push_back(n);
            n = ni(fwdpp::ts::nodes_preorder());
        }
    BOOST_REQUIRE_EQUAL(
        nodes == std::vector<fwdpp::ts::TS_NODE_INT>({ 6, 4, 0, 1, 5, 2, 3 }),
        true);
    auto nodes_from_fxn
        = fwdpp::ts::get_nodes(tv.tree(), fwdpp::ts::nodes_preorder());
    BOOST_REQUIRE_EQUAL(nodes == nodes_from_fxn, true);
}

BOOST_AUTO_TEST_SUITE_END()

