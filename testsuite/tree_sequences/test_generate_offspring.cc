#ifndef FWDPP_TEST_SUITE_MULTILOCUS_DETERMINISTIC_FIXTURE_HPP
#define FWDPP_TEST_SUITE_MULTILOCUS_DETERMINISTIC_FIXTURE_HPP

#include <iostream>
#include <config.h>
#include <cmath>
#include <queue>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/generate_offspring.hpp>
#include <fwdpp/ts/get_parent_ids.hpp>
#include <boost/test/unit_test.hpp>

#include "../fixtures/multilocus_fixture_deterministic.hpp"

BOOST_FIXTURE_TEST_CASE(check_multilocus_deterministic_fixture,
                        multilocus_fixture_deterministic)
{
    BOOST_REQUIRE_EQUAL(tables.genome_length(), nloci);
    BOOST_REQUIRE_EQUAL(pop.locus_boundaries.size(), nloci);
    BOOST_REQUIRE_EQUAL(mmodels.size(), nloci);
    BOOST_REQUIRE_EQUAL(intralocus_rec.size(), nloci);
    BOOST_REQUIRE_EQUAL(interlocus_rec.size(), nloci - 1);

    auto params = make_params();

    for (std::size_t i = 0; i < nloci; ++i)
        {
            auto k = params.generate_mutations[i](
                params.mutation_recycling_bin, pop.mutations);
            BOOST_REQUIRE_EQUAL(k.size(), 1);
            BOOST_REQUIRE_EQUAL(pop.mutations[k[0]].pos >= i, true);
            BOOST_REQUIRE_EQUAL(pop.mutations[k[0]].pos < i + 1, true);
            auto r = params.generate_breakpoints[i]();
            if (i % 2 == 0.)
                {
                    BOOST_REQUIRE_EQUAL(r.size(), 0);
                }
            else
                {
                    BOOST_REQUIRE_EQUAL(r.size(), 2);
                    BOOST_REQUIRE_EQUAL(r[0], i + 0.5);
                }
            if (i)
                {
                    auto ir = params.interlocus_recombination[i - 1]();
                    BOOST_REQUIRE_EQUAL(ir, i - 1);
                }
        }
}

BOOST_FIXTURE_TEST_CASE(test_transmission, multilocus_fixture_deterministic)
// This test applies mutations to a parent and tests transmission
// of those mutations into the offspring.
// This test is identical to one in unit/mlocusCrossoverTest.cc at face
// value.  The difference is that we are testing different machinery
// (fwdpp::ts::generate_offspring) and thus need a separate test
// to guard against future errors being introduced.
{
    mutate_parent();
    poptype::diploid_t offspring;
    auto params = make_params();
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);
    // Check transmission of mutations into offpring's FIRST gamete
    int locus = 0;
    bool expected_result = true;
    auto itr
        = std::find_if(begin(pop.gametes[offspring[locus].first].mutations),
                       end(pop.gametes[offspring[locus].first].mutations),
                       [this, locus](fwdpp::uint_t m) {
                           return pop.mutations[m].pos == locus;
                       });
    BOOST_REQUIRE_EQUAL(
        itr != end(pop.gametes[offspring[locus].first].mutations),
        expected_result);
    locus = 1;
    itr = std::find_if(begin(pop.gametes[offspring[locus].first].mutations),
                       end(pop.gametes[offspring[locus].first].mutations),
                       [this, locus](fwdpp::uint_t m) {
                           return pop.mutations[m].pos == locus;
                       });
    BOOST_REQUIRE_EQUAL(
        itr != end(pop.gametes[offspring[locus].first].mutations),
        expected_result);
    locus = 2;
    itr = std::find_if(begin(pop.gametes[offspring[locus].first].mutations),
                       end(pop.gametes[offspring[locus].first].mutations),
                       [this, locus](fwdpp::uint_t m) {
                           return pop.mutations[m].pos == locus;
                       });
    BOOST_REQUIRE_EQUAL(
        itr != end(pop.gametes[offspring[locus].first].mutations),
        expected_result);
    locus = 3;
    itr = std::find_if(begin(pop.gametes[offspring[locus].first].mutations),
                       end(pop.gametes[offspring[locus].first].mutations),
                       [this, locus](fwdpp::uint_t m) {
                           return pop.mutations[m].pos == locus;
                       });
    BOOST_REQUIRE_EQUAL(
        itr != end(pop.gametes[offspring[locus].first].mutations),
        expected_result);

    // Check transmission of mutations into offpring's SECOND gamete
    for (std::size_t i = 0; i < nloci; ++i)
        {
            locus = static_cast<int>(i);
            expected_result = false;
            itr = std::find_if(
                begin(pop.gametes[offspring[locus].second].mutations),
                end(pop.gametes[offspring[locus].second].mutations),
                [this, locus](fwdpp::uint_t m) {
                    return pop.mutations[m].pos == locus;
                });
            BOOST_REQUIRE_EQUAL(
                itr != end(pop.gametes[offspring[locus].second].mutations),
                expected_result);
        }
}

BOOST_FIXTURE_TEST_CASE(test_transmission_with_extra_variants,
                        multilocus_fixture_deterministic)
// This test applies mutations to a parent and tests transmission
// of those mutations into the offspring.
// This test is identical to one in unit/mlocusCrossoverTest.cc at face
// value.  The difference is that we are testing different machinery
// (fwdpp::ts::generate_offspring) and thus need a separate test
// to guard against future errors being introduced.
{
    mutate_parent2();
    poptype::diploid_t offspring;
    auto params = make_params();
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);
    for (std::size_t i = 0; i < nloci; ++i)
        {
            if (i % 2 == 0.)
                {
                    BOOST_REQUIRE_EQUAL(
                        pop.gametes[offspring[i].first].mutations.size(), 2);
                    BOOST_REQUIRE_EQUAL(
                        std::find_if(
                            begin(pop.gametes[offspring[i].first].mutations),
                            end(pop.gametes[offspring[i].first].mutations),
                            [this, i](fwdpp::uint_t k) {
                                return pop.mutations[k].pos == i;
                            })
                            != end(pop.gametes[offspring[i].first].mutations),
                        true);
                }
            else
                {
                    BOOST_REQUIRE_EQUAL(
                        pop.gametes[offspring[i].first].mutations.size(), 1);
                    BOOST_REQUIRE_EQUAL(
                        std::find_if(
                            begin(pop.gametes[offspring[i].first].mutations),
                            end(pop.gametes[offspring[i].first].mutations),
                            [this, i](fwdpp::uint_t k) {
                                return pop.mutations[k].pos == i;
                            })
                            != end(pop.gametes[offspring[i].first].mutations),
                        true);
                    BOOST_REQUIRE_EQUAL(
                        std::find_if(
                            begin(pop.gametes[offspring[i].first].mutations),
                            end(pop.gametes[offspring[i].first].mutations),
                            [this, i](fwdpp::uint_t k) {
                                return pop.mutations[k].pos == i + 0.51;
                            })
                            == end(pop.gametes[offspring[i].first].mutations),
                        true);
                }
        }
}

BOOST_FIXTURE_TEST_CASE(test_multilocus_determinisic_output,
                        multilocus_fixture_deterministic)
{
    auto params = make_params();
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);

    BOOST_CHECK_EQUAL(offspring.size(), nloci);
    BOOST_CHECK_EQUAL(pop.mutations.size(), 2 * nloci);
    // manually populate mcounts
    pop.mcounts.resize(pop.mutations.size());
    std::fill(begin(pop.mcounts), end(pop.mcounts), 0);
    for (std::size_t i = 0; i < nloci; ++i)
        {
            BOOST_CHECK_EQUAL(pop.gametes[offspring[i].first].mutations.size(),
                              0);
            BOOST_CHECK_EQUAL(
                pop.gametes[offspring[i].first].smutations.size(), 1);
            BOOST_CHECK_EQUAL(
                pop.gametes[offspring[i].second].mutations.size(), 0);
            BOOST_CHECK_EQUAL(
                pop.gametes[offspring[i].second].smutations.size(), 1);
            for (auto k : pop.gametes[offspring[i].first].smutations)
                {
                    pop.mcounts[k]++;
                }
            for (auto k : pop.gametes[offspring[i].second].smutations)
                {
                    pop.mcounts[k]++;
                }
            // Check variant positions
            for (auto k : pop.gametes[offspring[i].first].smutations)
                {
                    BOOST_CHECK_EQUAL(pop.mutations[k].pos >= i, true);
                    BOOST_CHECK_EQUAL(pop.mutations[k].pos < i + 1, true);
                }
        }
    BOOST_CHECK_EQUAL(std::accumulate(begin(pop.mcounts), end(pop.mcounts), 0),
                      2 * nloci);
    BOOST_CHECK_EQUAL(data_to_record.first.mutation_keys.size(), nloci);
    BOOST_CHECK_EQUAL(data_to_record.second.mutation_keys.size(), nloci);

    //TODO: test that the breakpoints contain the correct data!
}

BOOST_FIXTURE_TEST_CASE(test_multilocus_determinisic_output_swap_second,
                        multilocus_fixture_deterministic)
{
    auto params = make_params_swap_second();
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);
    BOOST_CHECK_EQUAL(offspring.size(), nloci);
    BOOST_CHECK_EQUAL(pop.mutations.size(), 2 * nloci);
    // manually populate mcounts
    pop.mcounts.resize(pop.mutations.size());
    std::fill(begin(pop.mcounts), end(pop.mcounts), 0);
    for (std::size_t i = 0; i < nloci; ++i)
        {
            BOOST_CHECK_EQUAL(pop.gametes[offspring[i].first].mutations.size(),
                              0);
            BOOST_CHECK_EQUAL(
                pop.gametes[offspring[i].first].smutations.size(), 1);
            BOOST_CHECK_EQUAL(
                pop.gametes[offspring[i].second].mutations.size(), 0);
            BOOST_CHECK_EQUAL(
                pop.gametes[offspring[i].second].smutations.size(), 1);
            for (auto k : pop.gametes[offspring[i].first].smutations)
                {
                    pop.mcounts[k]++;
                }
            for (auto k : pop.gametes[offspring[i].second].smutations)
                {
                    pop.mcounts[k]++;
                }
            // Check variant positions
            for (auto k : pop.gametes[offspring[i].first].smutations)
                {
                    BOOST_CHECK_EQUAL(pop.mutations[k].pos >= i, true);
                    BOOST_CHECK_EQUAL(pop.mutations[k].pos < i + 1, true);
                }
        }
    BOOST_CHECK_EQUAL(std::accumulate(begin(pop.mcounts), end(pop.mcounts), 0),
                      2 * nloci);
    BOOST_CHECK_EQUAL(data_to_record.first.mutation_keys.size(), nloci);
    BOOST_CHECK_EQUAL(data_to_record.second.mutation_keys.size(), nloci);

    std::vector<double> expected_breakpoints
        = { 1.5, 3.5, std::numeric_limits<double>::max() };
    for (auto e : expected_breakpoints)
        {
            auto itr = std::find(begin(data_to_record.first.breakpoints),
                                 end(data_to_record.first.breakpoints), e);
            auto itr2 = std::find(begin(data_to_record.second.breakpoints),
                                  end(data_to_record.second.breakpoints), e);
            BOOST_REQUIRE_EQUAL(itr != end(data_to_record.first.breakpoints),
                                true);
            BOOST_REQUIRE_EQUAL(itr2 != end(data_to_record.second.breakpoints),
                                true);
        }
}

BOOST_FIXTURE_TEST_CASE(test_multilocus_determinisic_table_recording,
                        multilocus_fixture_deterministic)
{
    auto params = make_params();
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);
    fwdpp::ts::TS_NODE_INT next_index = tables.node_table.size();
    auto p1d = fwdpp::ts::get_parent_ids(0, 0, data_to_record.first.swapped);
    auto p2d = fwdpp::ts::get_parent_ids(0, 1, data_to_record.second.swapped);
    tables.add_offspring_data(next_index++, data_to_record.first.breakpoints,
                              data_to_record.first.mutation_keys, p1d, 0, 1);
    tables.add_offspring_data(next_index++, data_to_record.second.breakpoints,
                              data_to_record.second.mutation_keys, p2d, 0, 1);
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size(), 8);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), 6);
}

BOOST_FIXTURE_TEST_CASE(
    test_multilocus_determinisic_table_recording_swap_second,
    multilocus_fixture_deterministic)
{
    auto params = make_params_swap_second();
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);
    fwdpp::ts::TS_NODE_INT next_index = tables.node_table.size();
    auto p1d = fwdpp::ts::get_parent_ids(0, 0, data_to_record.first.swapped);
    auto p2d = fwdpp::ts::get_parent_ids(0, 1, data_to_record.second.swapped);
    tables.add_offspring_data(next_index++, data_to_record.first.breakpoints,
                              data_to_record.first.mutation_keys, p1d, 0, 1);
    tables.add_offspring_data(next_index++, data_to_record.second.breakpoints,
                              data_to_record.second.mutation_keys, p2d, 0, 1);
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size(), 8);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), 6);
}

BOOST_FIXTURE_TEST_CASE(test_multilocus_determinisic_table_simplification,
                        multilocus_fixture_deterministic)
{
    auto params = make_params();
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params, offspring);
    fwdpp::ts::TS_NODE_INT next_index = tables.node_table.size();
    auto p1d = fwdpp::ts::get_parent_ids(0, 0, data_to_record.first.swapped);
    auto p2d = fwdpp::ts::get_parent_ids(0, 1, data_to_record.second.swapped);
    tables.add_offspring_data(next_index++, data_to_record.first.breakpoints,
                              data_to_record.first.mutation_keys, p1d, 0, 1);
    tables.add_offspring_data(next_index++, data_to_record.second.breakpoints,
                              data_to_record.second.mutation_keys, p2d, 0, 1);
    tables.sort_tables(pop.mutations);
    fwdpp::ts::table_simplifier simplifier(nloci);
    std::vector<fwdpp::ts::TS_NODE_INT> samples(
        { next_index - 2, next_index - 1 });
    auto rv = simplifier.simplify(tables, samples, pop.mutations);

    // All variants are new, so must survive simplification
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size(), 8);
}

#endif
