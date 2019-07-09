#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/marginal_tree_functions/children.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include "simple_table_collection.hpp"

int
naive_num_samples(const fwdpp::ts::marginal_tree &m, fwdpp::ts::TS_NODE_INT u)
{
    if (!m.advancing_sample_list())
        {
            throw std::runtime_error(
                "testsuite error: sample lists not being updated");
        }
    if (static_cast<std::size_t>(u) >= m.size())
        {
            throw std::runtime_error("testsuite error: node out of range");
        }
    int n = 0;
    auto l = m.left_sample[u];
    auto r = m.right_sample[u];
    while (true)
        {
            if (l == r)
                {
                    ++n;
                    break;
                }
            ++n;
            l = m.next_sample[l];
        }
    return n;
}

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
// NOTE: only works if all tips are samples.  This test can be used also
// as an independent count of the number of samples for this case.
{
    std::vector<fwdpp::ts::TS_NODE_INT> rv;
    get_tip(m, u, rv);
    std::sort(begin(rv), end(rv));
    return rv;
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
    tv = fwdpp::ts::tree_visitor(tables, samples,
                                 fwdpp::ts::update_samples_list(true));
    tv();
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
            tv = fwdpp::ts::tree_visitor(
                tables, samples, fwdpp::ts::update_samples_list(false));
            tv();
            // Throws b/c tv is not updating samples lists
            fwdpp::ts::num_samples(tv.tree(), 6);
        },
        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
