#include <fwdpp/ts/marginal_tree.hpp>
#include <boost/test/unit_test.hpp>
#include "simple_table_collection_infinite_sites.hpp"

BOOST_AUTO_TEST_SUITE()

BOOST_FIXTURE_TEST_CASE(test_construction_with_sample_groups,
                        simple_table_collection_infinite_sites)
{
    std::vector<fwdpp::ts::sample_group_map> groups;
    for (auto i : samples)
        {
            groups.emplace_back(i, 0);
        }
    fwdpp::ts::marginal_tree m(tables.node_table.size(), groups, true);
    BOOST_CHECK_EQUAL(m.sample_size(), groups.size());
    decltype(samples) scopy(m.samples_list_begin(), m.samples_list_end());
    BOOST_REQUIRE(samples == scopy);
}

BOOST_AUTO_TEST_SUITE_END()
