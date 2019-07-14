#include <iostream>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/marginal_tree_functions/roots.hpp>
#include <fwdpp/ts/marginal_tree_functions/nodes.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include <fwdpp/ts/generate_data_matrix.hpp>
#include "simple_table_collection_infinite_sites.hpp"

BOOST_AUTO_TEST_SUITE(test_generate_data_matrix)

BOOST_FIXTURE_TEST_CASE(infinite_sites, simple_table_collection_infinite_sites)
// There is exactly one mutation on each node, excluding the root.
{
    auto dm
        = fwdpp::ts::generate_data_matrix(tables, samples, true, true, false);
    BOOST_REQUIRE_EQUAL(dm.ncol, samples.size());
    auto rs = fwdpp::row_sums(dm);
    BOOST_REQUIRE(rs.second.empty());
    BOOST_REQUIRE_EQUAL(rs.first.size(), tables.node_table.size() - 1);
    auto roots = fwdpp::ts::get_roots(tv.tree());
    auto nodes = fwdpp::ts::get_nodes(tv.tree(), fwdpp::ts::nodes_preorder());
    std::vector<int> expected_sfs(samples.size(),
                                  0); //includes the fixed class;

    for (auto n : nodes)
        {
            // Exclude the root node from test b/c
            // the fixture doesn't put variants there.
            if (std::find(begin(roots), end(roots), n) == end(roots))
                {
                    auto lc = tv.tree().leaf_counts[n];
                    expected_sfs[lc - 1]++;
                }
        }
    for (std::size_t i = 0; i < expected_sfs.size(); ++i)
        {
            auto c = std::count(begin(rs.first), end(rs.first), i + 1);
            BOOST_CHECK_EQUAL(expected_sfs[i], c);
        }
}

BOOST_FIXTURE_TEST_CASE(simple_multiple_mutations_test,
                        simple_table_collection)
{
    tables.emplace_back_site(0.1, fwdpp::ts::default_ancestral_state);
    tables.emplace_back_mutation(4, 0lu, 0lu, fwdpp::ts::default_derived_state,
                                 true);
    auto d = fwdpp::ts::default_derived_state;
    d++;
    tables.emplace_back_mutation(4, 1lu, 0lu, d, true);
    auto dm
        = fwdpp::ts::generate_data_matrix(tables, samples, true, true, false);
    BOOST_REQUIRE_EQUAL(
        std::count(begin(dm.neutral.data), end(dm.neutral.data), d),
        fwdpp::ts::num_samples(tv.tree(), 4));
}

BOOST_AUTO_TEST_SUITE_END()
