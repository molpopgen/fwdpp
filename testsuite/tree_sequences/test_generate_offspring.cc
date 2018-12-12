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

/* Let's describe the fixture setup, expected outputs, etc.:
 *
 * The crossover functions in the fixture will assign the following numbers
 * of xover events w/in and b/w loci:
 *
 * Locus w/in b/w*
 * 0     0    N/A
 * 1     1    0
 * 2     0    1
 * 3     1    2
 *
 * *The entry for b/w loci refers to the number of xover events
 * b/w locus i-1 and i.
 *
 * An xover w/in locus i occurs at position i + 0.5
 *
 * Let's draw the 4 loci in two parent chromosomes:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++++++ ++++++++ ++++++++   
 * -------- -------- -------- --------   
 *
 * The way multi-locus recombination works is that we apply
 * each recombination event from left to right.  Internally, this
 * is not literally what happens, but it is our mental image of 
 * what we're trying to accomplish
 *
 * We will consider the case of NO INITIAL SWAPPING of either
 * parent gamete, and "build" the first offspring gamete.
 *
 * There are no x-overs w/in locus 1, so after processing it,
 * we should have:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++++++ ++++++++ ++++++++   
 *
 * Then, there is an xover w/in locus 1, giving:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- -------- --------   
 *
 * Then, 1 x-over event b/w loci 1 and 2:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++++++   
 *
 * No recombination events w/in locus 2:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++++++   
 *
 * Two recombination events b/w loci 2 and 3,
 * which means double x-over, which means 
 * no change in the system:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++++++   
 *
 * Finally, recombination w/in locus 3:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++----   
 *
 * Thus, the parent's first gamete gives 
 * the following interval's to the offspring:
 *
 * [0, 1.5)
 * [2, 3.5)
 *
 * The parent's second gamete passed on the following
 * intervals:
 *
 * [1.5, 2.0)
 * [3.5, 4.0)
 *
 * Thus, the breakpoints we need to send for tree 
 * sequence recording are:
 *
 * [1.5, 2.0, 3.5, DBL_MAX],
 *
 * These values are stored in "expected_breakpoints"
 * in the fixture.
 *
 * assuming that 4 is set as the "genome length"
 * of a table_collection.
 *
 * For the purpose of generating gametes, that "2.0"
 * breakpoint never needs to be stored, because parent 
 * gamete swapping, etc., takes care of that implicitly.
 *
 * Currently, we are FAILING to record it as a breakpoint
 * for ts recording, meaning that our tests below are 
 * WRONG.
 */

BOOST_FIXTURE_TEST_CASE(check_multilocus_deterministic_fixture,
                        multilocus_fixture_deterministic)
// This just sanity-checks various features of the fixture
{

    BOOST_REQUIRE_EQUAL(tables.genome_length(), nloci);
    BOOST_REQUIRE_EQUAL(pop.locus_boundaries.size(), nloci);
    BOOST_REQUIRE_EQUAL(mmodels.size(), nloci);
    BOOST_REQUIRE_EQUAL(intralocus_rec.size(), nloci);
    BOOST_REQUIRE_EQUAL(interlocus_rec.size(), nloci - 1);

    for (std::size_t i = 0; i < nloci; ++i)
        {
            auto k = params_no_swap.generate_mutations[i](
                params_no_swap.mutation_recycling_bin, pop.mutations);
            BOOST_REQUIRE_EQUAL(k.size(), 1);
            BOOST_REQUIRE_EQUAL(pop.mutations[k[0]].pos >= i, true);
            BOOST_REQUIRE_EQUAL(pop.mutations[k[0]].pos < i + 1, true);
            auto r = params_no_swap.generate_breakpoints[i]();
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
                    auto ir = params_no_swap.interlocus_recombination[i - 1]();
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
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap, offspring);
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
    BOOST_REQUIRE_NO_THROW(validate_mutations_positions_1(offspring));
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
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap, offspring);
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
    BOOST_REQUIRE_NO_THROW(validate_mutations_positions_2(offspring));
}

BOOST_FIXTURE_TEST_CASE(
    test_multilocus_determinisic_breakpoint_and_mutation_correctness,
    multilocus_fixture_deterministic)
{
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap, offspring);

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

BOOST_FIXTURE_TEST_CASE(
    test_multilocus_determinisic_output_breakponint_and_mutation_correctness_swap_second,
    multilocus_fixture_deterministic)
{
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_swap_second, offspring);
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
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap, offspring);
    fwdpp::ts::TS_NODE_INT next_index = tables.node_table.size();
    auto p1d = fwdpp::ts::get_parent_ids(0, 0, data_to_record.first.swapped);
    auto p2d = fwdpp::ts::get_parent_ids(0, 1, data_to_record.second.swapped);
    tables.add_offspring_data(next_index++, data_to_record.first.breakpoints,
                              data_to_record.first.mutation_keys, p1d, 0, 1);
    tables.add_offspring_data(next_index++, data_to_record.second.breakpoints,
                              data_to_record.second.mutation_keys, p2d, 0, 1);
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size(), 8);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), 8);
}

BOOST_FIXTURE_TEST_CASE(
    test_multilocus_determinisic_table_recording_swap_second,
    multilocus_fixture_deterministic)
{
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_swap_second, offspring);
    fwdpp::ts::TS_NODE_INT next_index = tables.node_table.size();
    auto p1d = fwdpp::ts::get_parent_ids(0, 0, data_to_record.first.swapped);
    auto p2d = fwdpp::ts::get_parent_ids(0, 1, data_to_record.second.swapped);
    tables.add_offspring_data(next_index++, data_to_record.first.breakpoints,
                              data_to_record.first.mutation_keys, p1d, 0, 1);
    tables.add_offspring_data(next_index++, data_to_record.second.breakpoints,
                              data_to_record.second.mutation_keys, p2d, 0, 1);
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size(), 8);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), 8);
}

BOOST_FIXTURE_TEST_CASE(test_multilocus_determinisic_table_simplification,
                        multilocus_fixture_deterministic)
{
    poptype::diploid_t offspring;
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap, offspring);
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
