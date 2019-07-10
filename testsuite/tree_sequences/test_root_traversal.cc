#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/decapitate.hpp>
#include <fwdpp/ts/marginal_tree_functions/roots.hpp>
#include "simple_table_collection.hpp"

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
	reset_visitor(false);
    auto r = fwdpp::ts::get_roots(tv.tree());
    BOOST_REQUIRE_EQUAL(r.size(), 2);
    BOOST_REQUIRE_EQUAL(r == expected_roots, true);

    fwdpp::ts::decapitate(tables, 1);
    expected_roots = { 2, 3, 4 };
	reset_visitor(false);
    r = fwdpp::ts::get_roots(tv.tree());
    std::sort(begin(r), end(r));
    BOOST_REQUIRE_EQUAL(r.size(), 3);
    BOOST_REQUIRE_EQUAL(r == expected_roots, true);

    fwdpp::ts::decapitate(tables, 2);
    expected_roots = { 0, 1, 2, 3 };
	reset_visitor(false);
    r = fwdpp::ts::get_roots(tv.tree());
    std::sort(begin(r), end(r));
    BOOST_REQUIRE_EQUAL(r.size(), 4);
    BOOST_REQUIRE_EQUAL(r == expected_roots, true);
}

BOOST_AUTO_TEST_SUITE_END()

