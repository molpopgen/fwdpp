#include <stdexcept>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/marginal_tree_functions/children.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include "simple_table_collection.hpp"

void
get_tip(const fwdpp::ts::marginal_tree &m, fwdpp::ts::TS_NODE_INT u,
        std::vector<fwdpp::ts::TS_NODE_INT> &samples)
{
    if (fwdpp::ts::num_children(m, u) > 0)
        {
            fwdpp::ts::process_children(
                m, u, true, [&m, &samples](fwdpp::ts::TS_NODE_INT x) {
                    get_tip(m, x, samples);
                });
            return;
        }
    samples.push_back(u);
}

std::vector<fwdpp::ts::TS_NODE_INT>
naive_get_samples(const fwdpp::ts::marginal_tree &m, fwdpp::ts::TS_NODE_INT u)
{
    if (!m.advancing_sample_list())
        {
            throw std::invalid_argument("sample lists are not being updated");
        }
    std::vector<fwdpp::ts::TS_NODE_INT> temp,
        samples_list(m.samples_list_begin(), m.samples_list_end()),
        intersection;
    std::sort(begin(samples_list), end(samples_list));
    get_tip(m, u, temp);
    std::sort(begin(temp), end(temp));
    std::set_intersection(begin(temp), end(temp), begin(samples_list),
                          end(samples_list), std::back_inserter(intersection));
    return intersection;
}

std::size_t
naive_num_samples(const fwdpp::ts::marginal_tree &m, fwdpp::ts::TS_NODE_INT u)
{
    auto s = naive_get_samples(m, u);
    return s.size();
}

BOOST_AUTO_TEST_SUITE(test_sample_traversal)

BOOST_FIXTURE_TEST_CASE(num_samples_simple_table_collection,
                        simple_table_collection)
{
    // NOTE: only works b/c no recombination, so node table size
    // is the number of nodes in the tree.  More generally,
    // we need to do a node traveral, but that hasn't been tested yet
    // (see below).
    for (std::size_t i = 0; i < tables.node_table.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(fwdpp::ts::num_samples(tv.tree(), i),
                                naive_num_samples(tv.tree(), i));
        }
}

BOOST_FIXTURE_TEST_CASE(test_sample_lists, simple_table_collection)
{
    for (std::size_t i = 0; i < tables.node_table.size(); ++i)
        {
            auto x = fwdpp::ts::get_samples(tv.tree(), i);
            auto y = naive_get_samples(tv.tree(), i);
            BOOST_CHECK_EQUAL(x.size(), y.size());
            BOOST_CHECK(fwdpp::ts::get_samples(tv.tree(), i)
                        == naive_get_samples(tv.tree(), i));
        }
}

BOOST_FIXTURE_TEST_CASE(test_subset_of_nodes_as_samples,
                        simple_table_collection)
{
    samples = { 0, 1, 3 };
    reset_visitor_and_samples(samples, true);
    BOOST_REQUIRE_EQUAL(naive_num_samples(tv.tree(), 6), samples.size());
    auto x = fwdpp::ts::get_samples(tv.tree(), 6);
    BOOST_CHECK(x.size() == samples.size());
    BOOST_CHECK(x == samples);
}

BOOST_FIXTURE_TEST_CASE(test_exception, simple_table_collection)
{
    // need a visitor that does NOT update samples lists
    BOOST_REQUIRE_THROW(
        {
            reset_visitor(false);
            // Throws b/c tv is not updating samples lists
            fwdpp::ts::num_samples(tv.tree(), 6);
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
