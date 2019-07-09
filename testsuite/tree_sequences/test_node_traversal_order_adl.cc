#include "preorder_adl.hpp"
#include "simple_table_collection.hpp"
#include <fwdpp/ts/marginal_tree_functions/nodes.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(test_adl)

BOOST_FIXTURE_TEST_CASE(test_preoder, simple_table_collection)
{
    auto n = fwdpp::ts::get_nodes(tv.tree(), nodes_preorder_test());
    auto n2 = fwdpp::ts::get_nodes(tv.tree(), fwdpp::ts::nodes_preorder());
    BOOST_REQUIRE(n == n2);
}

BOOST_AUTO_TEST_SUITE_END()

