#include <iostream>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection.hpp"
#include "simple_table_collection_polytomy.hpp"
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/decapitate.hpp>
#include <fwdpp/ts/marginal_tree_functions.hpp>

BOOST_AUTO_TEST_SUITE(test_traverse_children)

BOOST_FIXTURE_TEST_CASE(test_polytomy, simple_table_collection_polytomy)
{
    auto c = fwdpp::ts::get_children(tv.tree(), 7, true);
    decltype(c) expected = { 5, 6 };
    BOOST_REQUIRE(c == expected);

    c = fwdpp::ts::get_children(tv.tree(), 7, false);
    expected = { 6, 5 };
    BOOST_REQUIRE(c == expected);

    c = fwdpp::ts::get_children(tv.tree(), 5, false);
    expected = { 2, 1, 0 };
    BOOST_REQUIRE(c == expected);

    decltype(c) c2;
    fwdpp::ts::process_children(
        tv.tree(), 5, false,
        [&c2](fwdpp::ts::TS_NODE_INT u) { c2.push_back(u); });
    BOOST_REQUIRE(c == c2);
}

BOOST_FIXTURE_TEST_CASE(test_polytomy_decapitate,
                        simple_table_collection_polytomy)
{
    fwdpp::ts::decapitate(tables, 0.);
    tv = fwdpp::ts::tree_visitor(tables, samples,
                                 fwdpp::ts::update_samples_list(false));
    tv();
    auto c = fwdpp::ts::get_children(tv.tree(), 5, false);
    decltype(c) expected = { 2, 1, 0 };
    BOOST_REQUIRE(c == expected);
    BOOST_REQUIRE(fwdpp::ts::num_children(tv.tree(), 5) == 3);
    BOOST_REQUIRE(fwdpp::ts::num_children(tv.tree(), 0) == 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_root_traversal_simple_table_collection,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(test_roots)
{
    auto r = fwdpp::ts::get_roots(tv.tree());
    BOOST_REQUIRE_EQUAL(r.size(), 1);
    BOOST_REQUIRE_EQUAL(r[0], 6);
}

BOOST_AUTO_TEST_CASE(test_roots_after_decapitation)
{
    fwdpp::ts::decapitate(tables, 0.);
    std::vector<fwdpp::ts::TS_NODE_INT> expected_roots{ 4, 5 };
    tv = fwdpp::ts::tree_visitor(tables, samples,
                                 fwdpp::ts::update_samples_list(false));
    tv();
    auto r = fwdpp::ts::get_roots(tv.tree());
    BOOST_REQUIRE_EQUAL(r.size(), 2);
    BOOST_REQUIRE_EQUAL(r == expected_roots, true);

    fwdpp::ts::decapitate(tables, 1);
    expected_roots = { 2, 3, 4 };
    tv = fwdpp::ts::tree_visitor(tables, samples,
                                 fwdpp::ts::update_samples_list(false));
    tv();
    r = fwdpp::ts::get_roots(tv.tree());
    std::sort(begin(r), end(r));
    BOOST_REQUIRE_EQUAL(r.size(), 3);
    BOOST_REQUIRE_EQUAL(r == expected_roots, true);

    fwdpp::ts::decapitate(tables, 2);
    expected_roots = { 0, 1, 2, 3 };
    tv = fwdpp::ts::tree_visitor(tables, samples,
                                 fwdpp::ts::update_samples_list(false));
    tv();
    r = fwdpp::ts::get_roots(tv.tree());
    std::sort(begin(r), end(r));
    BOOST_REQUIRE_EQUAL(r.size(), 4);
    BOOST_REQUIRE_EQUAL(r == expected_roots, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_simple_table_collection_traversal,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_NO_THROW(
        { fwdpp::ts::node_iterator(tv.tree(), fwdpp::ts::nodes_preorder()); });
}

BOOST_AUTO_TEST_CASE(test_preorder_traversal)
{
    fwdpp::ts::node_iterator ni(tv.tree(), fwdpp::ts::nodes_preorder());
    std::vector<fwdpp::ts::TS_NODE_INT> nodes;
    auto n = ni();
    while (n != fwdpp::ts::TS_NULL_NODE)
        {
            nodes.push_back(n);
            n = ni();
        }
    BOOST_REQUIRE_EQUAL(
        nodes == std::vector<fwdpp::ts::TS_NODE_INT>({ 6, 4, 0, 1, 5, 2, 3 }),
        true);
    auto nodes_from_fxn
        = fwdpp::ts::get_nodes(tv.tree(), fwdpp::ts::nodes_preorder());
    BOOST_REQUIRE_EQUAL(nodes == nodes_from_fxn, true);
}

BOOST_AUTO_TEST_CASE(test_preorder_traversal_subtree)
{
    auto n = fwdpp::ts::get_nodes(tv.tree(), 5, fwdpp::ts::nodes_preorder());
    std::vector<fwdpp::ts::TS_NODE_INT> expected_nodes{ 5, 2, 3 };
    BOOST_REQUIRE_EQUAL(n.size(), expected_nodes.size());
    BOOST_REQUIRE_EQUAL(n == expected_nodes, true);
}

BOOST_AUTO_TEST_CASE(test_preorder_traversal_after_decaptitation)
{
    fwdpp::ts::decapitate(tables, 1.);
    tv = fwdpp::ts::tree_visitor(tables, samples,
                                 fwdpp::ts::update_samples_list(false));
    tv();
    auto nodes = fwdpp::ts::get_nodes(tv.tree(), fwdpp::ts::nodes_preorder());
    std::vector<fwdpp::ts::TS_NODE_INT> expected_nodes = { 4, 0, 1, 2, 3 };
    BOOST_REQUIRE_EQUAL(nodes.size(), expected_nodes.size());
    BOOST_REQUIRE_EQUAL(nodes.size() == expected_nodes.size(), true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_simple_table_collection_polytomy_traversal,
                         simple_table_collection_polytomy)

BOOST_AUTO_TEST_CASE(test_preorder_traversal)
{
    auto n = fwdpp::ts::get_nodes(tv.tree(), fwdpp::ts::nodes_preorder());
    BOOST_REQUIRE_EQUAL(n.size(), tables.num_nodes());
    std::vector<fwdpp::ts::TS_NODE_INT> expected_nodes{
        7, 5, 0, 1, 2, 6, 3, 4
    };
    BOOST_REQUIRE_EQUAL(n == expected_nodes, true);
}

BOOST_AUTO_TEST_CASE(test_preorder_traversal_subtree)
{
    auto n = fwdpp::ts::get_nodes(tv.tree(), 5, fwdpp::ts::nodes_preorder());
    std::vector<fwdpp::ts::TS_NODE_INT> expected_nodes{ 5, 0, 1, 2 };
    BOOST_REQUIRE_EQUAL(n == expected_nodes, true);
}

BOOST_AUTO_TEST_CASE(test_preorder_traversal_other_subtree)
{
    auto n = fwdpp::ts::get_nodes(tv.tree(), 6, fwdpp::ts::nodes_preorder());
    std::vector<fwdpp::ts::TS_NODE_INT> expected_nodes{ 6, 3, 4 };
    BOOST_REQUIRE_EQUAL(n == expected_nodes, true);
}

BOOST_AUTO_TEST_CASE(test_preorder_traversal_after_decaptitation)
{
    fwdpp::ts::decapitate(tables, 0.);
    tv = fwdpp::ts::tree_visitor(tables, samples,
                                 fwdpp::ts::update_samples_list(false));
    tv();
    auto nodes = fwdpp::ts::get_nodes(tv.tree(), fwdpp::ts::nodes_preorder());
    for (auto n : nodes)
        {
            std::cout << n << ' ';
        }
    std::cout << '\n';
    std::vector<fwdpp::ts::TS_NODE_INT> expected_nodes
        = { 5, 0, 1, 2, 6, 3, 4 };
    BOOST_REQUIRE_EQUAL(nodes.size(), expected_nodes.size());
    BOOST_REQUIRE_EQUAL(nodes.size() == expected_nodes.size(), true);
}

BOOST_AUTO_TEST_SUITE_END()

