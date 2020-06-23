#include <fwdpp/ts/tree_visitor.hpp>
#include <boost/test/unit_test.hpp>
#include "empty_table_collection.hpp"
#include "simple_table_collection.hpp"
#include "simple_table_collection_polytomy.hpp"

#include "wfevolve_table_collection.hpp"
#include "tskit_utils.hpp"

namespace
{
    struct ll_wf_fixture
    {
        unsigned N;
        unsigned nsteps;
        fwdpp::ts::std_table_collection tables;
        wf_simulation_results results;
        std::vector<fwdpp::ts::sample_group_map> samples;

        std::vector<fwdpp::ts::sample_group_map>
        fill_samples(const wf_simulation_results& results)
        {
            std::vector<fwdpp::ts::sample_group_map> rv;
            for (auto& p : results.alive_individuals)
                {
                    for (auto n : p.nodes)
                        {
                            if (n < static_cast<fwdpp::ts::table_index_t>(
                                    results.alive_individuals.size()))
                                {
                                    rv.emplace_back(n, 0);
                                }
                            else
                                {
                                    rv.emplace_back(n, 1);
                                }
                        }
                }
            return rv;
        }

        explicit ll_wf_fixture(wf_simulation_params p)
            : N{p.N}, nsteps{p.nsteps}, tables{1.},
              results{wfevolve_table_collection(p.seed, p.N, p.nsteps, p.psurvival,
                                                p.rho, p.simplification_interval, true,
                                                true, false, empty_policies{}, tables)},
              samples{fill_samples(results)}
        {
            tables.build_indexes();
        }
    };

    bool
    validate_sample_group_counts(const std::vector<fwdpp::ts::sample_group_map>& samples,
                                 int group, int observed_number)
    {
        int expected_number = 0;
        for (auto& s : samples)
            {
                if (s.group == group)
                    {
                        ++expected_number;
                    }
            }
        return expected_number == observed_number;
    }

    struct wf_fixture_no_rec_non_overlapping
    {
        ll_wf_fixture ll;
        wf_fixture_no_rec_non_overlapping() : ll({42, 100, 1000, 0., 0., 100})
        {
        }
    };

    struct wf_fixture_with_rec_non_overlapping
    {
        ll_wf_fixture ll;
        wf_fixture_with_rec_non_overlapping() : ll({42, 100, 1000, 0., 100., 100})
        {
        }
    };

}

BOOST_FIXTURE_TEST_SUITE(test_tree_visitor_with_simple_table_collection,
                         simple_table_collection)

BOOST_AUTO_TEST_CASE(test_unindexed_tables)
{
    tables.input_left.clear();
    tables.output_right.clear();
    BOOST_REQUIRE_THROW(fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection> tv(
                            tables, samples, fwdpp::ts::update_samples_list(0)),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_tree_visitor_with_empty_table_collection,
                         empty_table_collection)

BOOST_AUTO_TEST_CASE(test_construction)
{
    BOOST_REQUIRE_THROW(fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection> tv(
                            tables, empty_samples, fwdpp::ts::update_samples_list(0)),
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
    tv = fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>(
        tables, groups, fwdpp::ts::update_samples_list(1));
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

BOOST_FIXTURE_TEST_CASE(test_group_labels_from_simulation,
                        wf_fixture_no_rec_non_overlapping)
{
    auto tv = fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>(
        ll.tables, ll.samples, fwdpp::ts::update_samples_list(1));
    std::vector<int> nodes(ll.tables.num_nodes());
    std::iota(begin(nodes), end(nodes), 0);
    while (tv())
        {
            std::vector<int> is_sample(nodes.size(), 0);
            int n0 = 0, n1 = 0;
            for (auto& p : ll.results.alive_individuals)
                {
                    for (auto n : p.nodes)
                        {
                            if (n < static_cast<fwdpp::ts::table_index_t>(
                                    ll.results.alive_individuals.size()))
                                {
                                    BOOST_REQUIRE_EQUAL(tv.tree().sample_group(n), 0);
                                    is_sample[n] = 1;
                                    ++n0;
                                }
                            else
                                {
                                    BOOST_REQUIRE_EQUAL(tv.tree().sample_group(n), 1);
                                    is_sample[n] = 1;
                                    ++n0;
                                }
                        }
                }
            validate_sample_group_counts(ll.samples, 0, n0);
            validate_sample_group_counts(ll.samples, 1, n1);
            for (std::size_t i = 0; i < is_sample.size(); ++i)
                {
                    if (is_sample[i] == 0)
                        {
                            BOOST_REQUIRE(tv.tree().sample_group(i) < 0);
                        }
                }
        }
}

BOOST_FIXTURE_TEST_CASE(test_group_labels_from_simulation_with_recombination,
                        wf_fixture_with_rec_non_overlapping)
{
    auto tv = fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>(
        ll.tables, ll.samples, fwdpp::ts::update_samples_list(1));
    std::vector<int> nodes(ll.tables.num_nodes());
    std::iota(begin(nodes), end(nodes), 0);
    while (tv())
        {
            std::vector<int> is_sample(nodes.size(), 0);
            int n0 = 0, n1 = 0;
            for (auto& p : ll.results.alive_individuals)
                {
                    for (auto n : p.nodes)
                        {
                            if (n < static_cast<fwdpp::ts::table_index_t>(
                                    ll.results.alive_individuals.size()))
                                {
                                    BOOST_REQUIRE_EQUAL(tv.tree().sample_group(n), 0);
                                    is_sample[n] = 1;
                                    ++n0;
                                }
                            else
                                {
                                    BOOST_REQUIRE_EQUAL(tv.tree().sample_group(n), 1);
                                    is_sample[n] = 1;
                                    ++n1;
                                }
                        }
                }
            validate_sample_group_counts(ll.samples, 0, n0);
            validate_sample_group_counts(ll.samples, 1, n1);
            for (std::size_t i = 0; i < is_sample.size(); ++i)
                {
                    if (is_sample[i] == 0)
                        {
                            BOOST_REQUIRE(tv.tree().sample_group(i) < 0);
                        }
                }
        }
}

BOOST_AUTO_TEST_SUITE_END()
