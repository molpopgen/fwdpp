#include <fwdpp/ts/marginal_tree.hpp>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection_infinite_sites.hpp"

BOOST_AUTO_TEST_SUITE()

BOOST_FIXTURE_TEST_CASE(test_construction_with_sample_groups,
                        simple_table_collection_infinite_sites)
{
    std::vector<fwdpp::ts::sample_group_map<fwdpp::ts::table_collection::id_type>>
        groups;
    for (auto i : samples)
        {
            groups.emplace_back(i, 0);
        }
    fwdpp::ts::marginal_tree<fwdpp::ts::table_collection::id_type> m(tables.nodes.size(),
                                                                     groups, true);
    BOOST_CHECK_EQUAL(m.sample_size(), groups.size());
    decltype(samples) scopy(m.samples_list_begin(), m.samples_list_end());
    BOOST_REQUIRE(samples == scopy);
}

BOOST_FIXTURE_TEST_CASE(test_construction_with_sample_groups_and_preserved_nodes,
                        simple_table_collection_infinite_sites)
{
    std::vector<fwdpp::ts::sample_group_map<fwdpp::ts::table_collection::id_type>>
        groups;
    for (auto i : samples)
        {
            groups.emplace_back(i, 0);
        }
    std::vector<fwdpp::ts::table_collection::id_type> preserved_nodes;
    preserved_nodes.push_back(
        static_cast<fwdpp::ts::table_collection::id_type>(samples.size()));
    fwdpp::ts::marginal_tree<fwdpp::ts::table_collection::id_type> m(
        tables.nodes.size(), groups, preserved_nodes, true);
    BOOST_CHECK_EQUAL(m.sample_size(), groups.size() + preserved_nodes.size());
    decltype(samples) scopy(m.samples_list_begin(), m.samples_list_end());
    samples.insert(end(samples), begin(preserved_nodes), end(preserved_nodes));
    BOOST_REQUIRE(samples == scopy);
}

BOOST_FIXTURE_TEST_CASE(test_sample_preserved_node_overlap,
                        simple_table_collection_infinite_sites)
{
    BOOST_REQUIRE_THROW(
        {
            fwdpp::ts::marginal_tree<fwdpp::ts::table_collection::id_type> m(
                tables.nodes.size(), samples, samples, true);
        },
        fwdpp::ts::samples_error);
}
BOOST_AUTO_TEST_SUITE_END()
