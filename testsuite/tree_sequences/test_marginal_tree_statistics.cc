#include <iostream>
#include <fwdpp/ts/marginal_tree_functions/statistics.hpp>
#include <fwdpp/ts/decapitate.hpp>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection.hpp"
#include "simple_table_collection_polytomy.hpp"

double
naive_branch_length(const fwdpp::ts::marginal_tree &m,
                    const std::vector<fwdpp::ts::TS_NODE_INT> &samples,
                    const std::vector<fwdpp::ts::node> &nodes,
                    bool scale_by_length)
// Sample-to-root implementation requiring O(nnodes) extra memory,
// and makes multiple passes through same ancestral nodes.
{
    double ttime = 0;
    std::vector<std::int8_t> processed(nodes.size(), 0);
    for (auto s : samples)
        {
            auto u = s;
            while (u != fwdpp::ts::TS_NULL_NODE)
                {
                    auto p = m.parents[u];
                    if (!processed[u] && p != fwdpp::ts::TS_NULL_NODE)
                        {
                            ttime += (nodes[u].time - nodes[p].time);
                        }
                    processed[u] = 1;
                    u = p;
                }
        }
    if (scale_by_length)
        {
            ttime *= (m.right - m.left);
        }
    return ttime;
}

BOOST_AUTO_TEST_SUITE(test_marginal_tree_statistics)

BOOST_FIXTURE_TEST_CASE(test_total_time_simple, simple_table_collection)
{
    auto ttime = fwdpp::ts::total_time(tv.tree(), tables.node_table, true);
    BOOST_REQUIRE_EQUAL(ttime, total_time);
    BOOST_REQUIRE_EQUAL(ttime, naive_branch_length(tv.tree(), samples,
                                                   tables.node_table, true));
}

BOOST_FIXTURE_TEST_CASE(test_total_time_polytomy,
                        simple_table_collection_polytomy)
{
    auto ttime = fwdpp::ts::total_time(tv.tree(), tables.node_table, true);
    BOOST_REQUIRE_EQUAL(ttime, total_time);
}

BOOST_FIXTURE_TEST_CASE(test_total_time_polytomy_decapitate,
                        simple_table_collection_polytomy)
{
    fwdpp::ts::decapitate(tables, 0.);
    tv = fwdpp::ts::tree_visitor(tables, samples);
    tv(std::false_type(), std::false_type());
    auto ttime = fwdpp::ts::total_time(tv.tree(), tables.node_table, true);
    BOOST_REQUIRE(ttime != total_time);
    BOOST_REQUIRE_EQUAL(ttime, naive_branch_length(tv.tree(), samples,
                                                   tables.node_table, true));
}

BOOST_AUTO_TEST_SUITE_END()

