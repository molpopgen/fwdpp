#include <iostream>
#include <boost/test/unit_test.hpp>
#include <fwdpp/ts/std_table_collection.hpp>
#include "wfevolve_table_collection.hpp"
#include "tskit_utils.hpp"

namespace
{
    struct ll_wf_fixture
    {
        unsigned N;
        unsigned nsteps;
        fwdpp::ts::std_table_collection tables_buffering, tables_sorting;
        wf_simulation_results results_buffering, results_sorting;
        std::vector<int> samples_buffering, samples_sorting;

        std::vector<int>
        fill_samples(const fwdpp::ts::std_table_collection& tables,
                     const wf_simulation_results& results)
        {
            std::vector<int> rv(tables.num_nodes(), 0);
            const auto update = [](auto& v, auto n) {
                if (n < 0)
                    {
                        throw std::runtime_error("node < 0");
                    }
                if (n >= static_cast<decltype(n)>(v.size()))
                    {
                        throw std::runtime_error("node out of bounds");
                    }
                v[n] = 1;
            };

            for (auto& p : results.alive_individuals)
                {
                    for (auto n : p.nodes)
                        {
                            update(rv, n);
                        }
                }
            for (auto n : results.preserved_nodes)
                {
                    update(rv, n);
                }
            return rv;
        }

        explicit ll_wf_fixture(wf_simulation_params p)
            : N{p.N}, nsteps{p.nsteps}, tables_buffering{1.}, tables_sorting{1.},
              results_buffering{wfevolve_table_collection(
                  p.seed, p.N, p.nsteps, p.psurvival, p.rho, p.simplification_interval,
                  true, true, false, empty_policies{}, tables_buffering)},
              results_sorting{wfevolve_table_collection(
                  p.seed, p.N, p.nsteps, p.psurvival, p.rho, p.simplification_interval,
                  false, false, false, empty_policies{}, tables_sorting)},
              samples_buffering{fill_samples(tables_buffering, results_buffering)},
              samples_sorting{fill_samples(tables_sorting, results_sorting)}
        {
        }
    };

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
        wf_fixture_with_rec_non_overlapping() : ll({4298764, 100, 3000, 0., 5., 100})
        {
        }
    };
}

BOOST_FIXTURE_TEST_SUITE(no_recombination_non_overlapping,
                         wf_fixture_no_rec_non_overlapping)

BOOST_AUTO_TEST_CASE(test_edge_table_size)
{
    BOOST_REQUIRE_EQUAL(ll.tables_buffering.num_nodes(), ll.tables_sorting.num_nodes());
}

BOOST_AUTO_TEST_CASE(test_node_table_size)
{
    BOOST_REQUIRE_EQUAL(ll.tables_buffering.num_edges(), ll.tables_sorting.num_edges());
}

BOOST_AUTO_TEST_CASE(test_node_table_time_equality)
{
    for (std::size_t i = 0; i < ll.tables_buffering.num_nodes(); ++i)
        {
            BOOST_REQUIRE_EQUAL(ll.tables_buffering.nodes[i].time,
                                ll.tables_sorting.nodes[i].time);
        }
}

BOOST_AUTO_TEST_CASE(test_kc_distance)
{
    auto tsk_tables_buffer = dump_table_collection_to_tskit(
        ll.tables_buffering, ll.nsteps, ll.samples_buffering);
    auto tsk_tables_sort = dump_table_collection_to_tskit(ll.tables_buffering, ll.nsteps,
                                                          ll.samples_sorting);
    tsk_treeseq_wrapper treeseq_buffer(tsk_tables_buffer.get());
    tsk_treeseq_wrapper treeseq_sort(tsk_tables_buffer.get());
    double kc_distance = std::numeric_limits<double>::quiet_NaN();
    BOOST_REQUIRE_NO_THROW({
        auto rv = tsk_treeseq_kc_distance(treeseq_buffer.get(), treeseq_sort.get(), 0.,
                                          &kc_distance);
        handle_tskit_return_code(rv);
    });
    BOOST_REQUIRE_EQUAL(kc_distance, 0.);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(with_recombination_non_overlapping,
                         wf_fixture_with_rec_non_overlapping)

BOOST_AUTO_TEST_CASE(test_edge_table_size)
{
    BOOST_REQUIRE_EQUAL(ll.tables_buffering.num_nodes(), ll.tables_sorting.num_nodes());
}

BOOST_AUTO_TEST_CASE(test_node_table_size)
{
    BOOST_REQUIRE_EQUAL(ll.tables_buffering.num_edges(), ll.tables_sorting.num_edges());
}

BOOST_AUTO_TEST_CASE(test_node_table_time_equality)
{
    for (std::size_t i = 0; i < ll.tables_buffering.num_nodes(); ++i)
        {
            BOOST_REQUIRE_EQUAL(ll.tables_buffering.nodes[i].time,
                                ll.tables_sorting.nodes[i].time);
        }
}

BOOST_AUTO_TEST_CASE(test_kc_distance)
{
    auto tsk_tables_buffer = dump_table_collection_to_tskit(
        ll.tables_buffering, ll.nsteps, ll.samples_buffering);
    auto tsk_tables_sort = dump_table_collection_to_tskit(ll.tables_buffering, ll.nsteps,
                                                          ll.samples_sorting);
    tsk_treeseq_wrapper treeseq_buffer(tsk_tables_buffer.get());
    tsk_treeseq_wrapper treeseq_sort(tsk_tables_buffer.get());
    BOOST_REQUIRE(treeseq_buffer.treeseq.num_trees > 0);
    BOOST_REQUIRE_EQUAL(treeseq_buffer.treeseq.num_trees,
                        treeseq_sort.treeseq.num_trees);
    double kc_distance = std::numeric_limits<double>::quiet_NaN();
    BOOST_REQUIRE_NO_THROW({
        auto rv = tsk_treeseq_kc_distance(treeseq_buffer.get(), treeseq_sort.get(), 0.,
                                          &kc_distance);
        handle_tskit_return_code(rv);
    });
    BOOST_REQUIRE_EQUAL(kc_distance, 0.);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(overlapping)

BOOST_AUTO_TEST_CASE(test_no_recombination_overlapping)
{
    fwdpp::ts::std_table_collection tables_buffering{1.};
    auto tables_sorting{tables_buffering};

    unsigned seed = 42;
    unsigned N = 1000;
    unsigned nsteps = 5 * N;
    double psurvival = 0.5;
    double rho = 0.;
    unsigned simplification_interval = 100;

    BOOST_REQUIRE_NO_THROW({
        wfevolve_table_collection(seed, N, nsteps, psurvival, rho,
                                  simplification_interval, true, true, false,
                                  empty_policies{}, tables_buffering);
    });
    BOOST_REQUIRE_NO_THROW({
        wfevolve_table_collection(seed, N, nsteps, psurvival, rho,
                                  simplification_interval, false, false, false,
                                  empty_policies{}, tables_sorting);
    });
    BOOST_REQUIRE_EQUAL(tables_buffering.num_nodes(), tables_sorting.num_nodes());
    BOOST_REQUIRE_EQUAL(tables_buffering.num_edges(), tables_sorting.num_edges());
    for (std::size_t i = 0; i < tables_buffering.num_nodes(); ++i)
        {
            BOOST_REQUIRE_EQUAL(tables_buffering.nodes[i].time,
                                tables_sorting.nodes[i].time);
        }
    // NOTE: KC-distance bits in tskit don't support overlapping
    // generations at the moment.  You get a "unsimplified tree
    // with unary nodes error"

    //auto tsk_tables_buffer = dump_table_collection_to_tskit(tables_buffering, nsteps, N);
    //auto tsk_tables_sort = dump_table_collection_to_tskit(tables_sorting, nsteps, N);
    //tsk_treeseq_wrapper treeseq_buffer(tsk_tables_buffer.get());
    //tsk_treeseq_wrapper treeseq_sort(tsk_tables_buffer.get());
    //double kc_distance = std::numeric_limits<double>::quiet_NaN();
    //BOOST_REQUIRE_NO_THROW({
    //    auto rv = tsk_treeseq_kc_distance(treeseq_buffer.get(), treeseq_sort.get(), 0.,
    //                                      &kc_distance);
    //    handle_tskit_return_code(rv);
    //});
    //BOOST_REQUIRE_EQUAL(kc_distance, 0.);
}

BOOST_AUTO_TEST_SUITE_END()

