#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/decapitate.hpp>
#include <fwdpp/ts/marginal_tree_functions/children.hpp>
#include "simple_table_collection.hpp"
#include "simple_table_collection_polytomy.hpp"

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
        [&c2](fwdpp::ts::table_index_t u) { c2.push_back(u); });
    BOOST_REQUIRE(c == c2);
}

BOOST_FIXTURE_TEST_CASE(test_polytomy_decapitate,
                        simple_table_collection_polytomy)
{
    fwdpp::ts::decapitate(tables, 0., false);
	reset_visitor(false);
    auto c = fwdpp::ts::get_children(tv.tree(), 5, false);
    decltype(c) expected = { 2, 1, 0 };
    BOOST_REQUIRE(c == expected);
    BOOST_REQUIRE(fwdpp::ts::num_children(tv.tree(), 5) == 3);
    BOOST_REQUIRE(fwdpp::ts::num_children(tv.tree(), 0) == 0);
}

BOOST_AUTO_TEST_SUITE_END()

