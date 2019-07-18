#include <fwdpp/ts/tree_visitor.hpp>
#include <boost/test/unit_test.hpp>
#include "empty_table_collection.hpp"
#include "simple_table_collection.hpp"
#include "simple_table_collection_polytomy.hpp"

BOOST_FIXTURE_TEST_SUITE(test_tree_visitor_with_simple_table_collection,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(test_construction)
{
    // NOTE: this is a "cheat", as we cannot get this far
    // if this doesn't work, as the fixture constructor
    // does this.
    BOOST_REQUIRE_NO_THROW({
        fwdpp::ts::tree_visitor(tables, samples,
                                fwdpp::ts::update_samples_list(0));
    });
}

BOOST_AUTO_TEST_CASE(test_unindexed_tables)
{
    tables.input_left.clear();
    tables.output_right.clear();
    BOOST_REQUIRE_THROW(
        fwdpp::ts::tree_visitor tv(tables, samples,
                                   fwdpp::ts::update_samples_list(0)),
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_tree_visitor_with_empty_table_collection,
                         empty_table_collection)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_THROW(
        fwdpp::ts::tree_visitor tv(tables, empty_samples,
                                   fwdpp::ts::update_samples_list(0)),
        fwdpp::ts::samples_error);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_sample_groups)

BOOST_FIXTURE_TEST_CASE(test_group_labels, simple_table_collection_polytomy)
{
    std::vector<fwdpp::ts::sample_group_map> groups;
    for (auto i : samples)
        {
            groups.emplace_back(i, 0);
        }
    groups[2].group = 1;
    tv = fwdpp::ts::tree_visitor(tables, groups,
                                 fwdpp::ts::update_samples_list(1));
    tv();
    auto c = tv.tree().left_sample[5];
    BOOST_REQUIRE(c == 0);
    BOOST_REQUIRE(tv.tree().sample_group(c) == 0);
    c = tv.tree().next_sample[c];
    BOOST_REQUIRE(c == 1);
    BOOST_REQUIRE(tv.tree().sample_group(c) == 0);
    c = tv.tree().next_sample[c];
    BOOST_REQUIRE(c == 2);
    BOOST_REQUIRE(tv.tree().sample_group(c) == 1);
}

BOOST_AUTO_TEST_SUITE_END()
