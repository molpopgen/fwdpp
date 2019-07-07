#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/decapitate.hpp>
#include "simple_table_collection.hpp"

BOOST_FIXTURE_TEST_SUITE(test_decaptiate_table_collection,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(remove_root)
{
    auto n = tables.num_nodes();
    auto e = tables.edge_table.size();
    fwdpp::ts::decapitate(tables, 0.0);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 1);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), e - 2);
}

BOOST_AUTO_TEST_CASE(remove_root_and_node_5)
{
    auto n = tables.num_nodes();
    auto e = tables.edge_table.size();
    fwdpp::ts::decapitate(tables, 1.0);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 2);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), e - 4);
}

BOOST_AUTO_TEST_CASE(remove_root_thru_node_4)
{
    auto n = tables.num_nodes();
    auto e = tables.edge_table.size();
    fwdpp::ts::decapitate(tables, 2.0);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 3);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), e - 6);
}

BOOST_AUTO_TEST_CASE(remove_all_nodes)
{
    fwdpp::ts::decapitate(tables, 4.0);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), 0);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
