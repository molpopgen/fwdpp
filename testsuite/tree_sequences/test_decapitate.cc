#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/decapitate.hpp>
#include <fwdpp/ts/site.hpp>
#include <fwdpp/ts/mutation_record.hpp>
#include "simple_table_collection.hpp"

BOOST_FIXTURE_TEST_SUITE(test_decaptiate_table_collection,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(remove_root)
{
    auto n = tables.num_nodes();
    auto e = tables.edges.size();
    fwdpp::ts::decapitate(tables, 0.0, false);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 1);
    BOOST_REQUIRE_EQUAL(tables.edges.size(), e - 2);
}

BOOST_AUTO_TEST_CASE(remove_root_and_node_5)
{
    auto n = tables.num_nodes();
    auto e = tables.edges.size();
    fwdpp::ts::decapitate(tables, 1.0, false);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 2);
    BOOST_REQUIRE_EQUAL(tables.edges.size(), e - 4);
}

BOOST_AUTO_TEST_CASE(remove_root_and_node_5_with_mutation)
{
    tables.sites.emplace_back(fwdpp::ts::site{0.25, 0});
    tables.sites.emplace_back(fwdpp::ts::site{0.666, 0});
    tables.mutations.emplace_back(fwdpp::ts::mutation_record{4, 0, 0, 1, true});
    tables.mutations.emplace_back(fwdpp::ts::mutation_record{5, 1, 0, 1, true});
    auto n = tables.num_nodes();
    auto e = tables.edges.size();
    fwdpp::ts::decapitate(tables, 1.0, true);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 2);
    BOOST_REQUIRE_EQUAL(tables.edges.size(), e - 4);
    BOOST_REQUIRE_EQUAL(tables.mutations.size(), 1);
    BOOST_REQUIRE_EQUAL(tables.sites[0].position, 0.25);
    BOOST_REQUIRE_EQUAL(tables.mutations[0].site, 0);
    BOOST_REQUIRE_EQUAL(tables.sites.size(), 1);
}

BOOST_AUTO_TEST_CASE(remove_root_and_node_5_with_mutation_dont_prune_mutation)
{
    tables.sites.emplace_back(fwdpp::ts::site{0.666, 0});
    tables.mutations.emplace_back(fwdpp::ts::mutation_record{5, 0, 0, 1, true});
    auto n = tables.num_nodes();
    auto e = tables.edges.size();
    fwdpp::ts::decapitate(tables, 1.0, false);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 2);
    BOOST_REQUIRE_EQUAL(tables.edges.size(), e - 4);
    BOOST_REQUIRE_EQUAL(tables.mutations.size(), 1);
    BOOST_REQUIRE_EQUAL(tables.sites.size(), 1);
}

BOOST_AUTO_TEST_CASE(remove_root_thru_node_4)
{
    auto n = tables.num_nodes();
    auto e = tables.edges.size();
    fwdpp::ts::decapitate(tables, 2.0, false);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), n - 3);
    BOOST_REQUIRE_EQUAL(tables.edges.size(), e - 6);
}

BOOST_AUTO_TEST_CASE(remove_all_nodes)
{
    fwdpp::ts::decapitate(tables, 4.0, false);
    BOOST_REQUIRE_EQUAL(tables.num_nodes(), 0);
    BOOST_REQUIRE_EQUAL(tables.edges.size(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
