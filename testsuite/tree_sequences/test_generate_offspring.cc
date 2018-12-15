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
#include <fwdpp/ts/count_mutations.hpp>
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
 * assuming that 4 is set as the "genome length"
 * of a table_collection.
 *
 * These values are stored in "expected_breakpoints"
 * in the fixture.
 *
 *
 * For the purpose of generating gametes, that "2.0"
 * breakpoint never needs to be stored, because parent 
 * gamete swapping, etc., takes care of that implicitly.
 *
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
    BOOST_REQUIRE_EQUAL(data_to_record.first.swapped, 0);
    BOOST_REQUIRE_EQUAL(data_to_record.second.swapped, 1);
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

//NOTE: test are believed to be correct to this point.

BOOST_FIXTURE_TEST_CASE(test_multilocus_determinisic_table_simplification,
                        multilocus_fixture_deterministic)
//Generate a single offspring after mutating parents and then test the
//transmission of mutations through simplification w.r.to that single
//offspring
{
    // Calling mutate_parent2 places mutations at positions i and i+0.51
    // for i in [0,4) on all "first" gametes of diploid 0, so that is two mutations per locus.
    mutate_parent2();
    // Using params_no swap means that the following mutations will be passed from node zero
    // to node 2N:
    // 0, 0.51, 1, 2, 2.51, 3.0
    // Mutations at the following positions must be simplified out of the mutation table:
    // 1.51, 2.51
    // This, we have the following expectations:
    std::vector<double> muts_on_node_0
        = { 0., 0.51, 1., 2., 2.51, 3. }, //node_0 means post-simplification!
        muts_lost = { 1.51, 3.51 };
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size(), 8);
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

    // We have simplified to the two extant samples, which are the
    // to genomes of "offspring". Thus, there are zero transmission
    // events required to describe the ancestry, so the edge table is
    // empty.
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), 0);
    BOOST_REQUIRE_EQUAL(tables.node_table.size(), 2);

    for (auto& n : tables.node_table)
        {
            BOOST_REQUIRE_EQUAL(n.time, 1.);
        }

    // Check the status of mutations transmitted from parents to offspring
    for (auto p : muts_on_node_0)
        {
            bool found = false;
            for (auto m : tables.mutation_table)
                {
                    if (m.node == 0 && pop.mutations[m.key].pos == p)
                        {
                            found = true;
                        }
                }
            BOOST_REQUIRE_EQUAL(found, true);
        }

    for (auto p : muts_lost)
        {
            bool found = false;
            for (auto m : tables.mutation_table)
                {
                    if (pop.mutations[m.key].pos == p)
                        {
                            found = true;
                        }
                }
            BOOST_REQUIRE_EQUAL(found, false);
        }
    // Count numbers of new mutations
    int new_node_zero = 0, new_node_one = 0;
    for (auto& m : tables.mutation_table)
        {
            if (m.node == 0 && pop.mutations[m.key].g == 1)
                {
                    ++new_node_zero;
                }
            if (m.node == 1 && pop.mutations[m.key].g == 1)
                {
                    ++new_node_one;
                }
        }
    BOOST_REQUIRE_EQUAL(new_node_zero, 4);
    BOOST_REQUIRE_EQUAL(new_node_one, 4);
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size() - new_node_one
                            - new_node_zero,
                        muts_on_node_0.size());
}

BOOST_FIXTURE_TEST_CASE(
    test_multilocus_determinisic_table_simplification_then_count_mutations,
    multilocus_fixture_deterministic)
//Same as previous test, then we count mutations
{
    // Calling mutate_parent2 places mutations at positions i and i+0.51
    // for i in [0,4) on all "first" gametes of diploid 0, so that is two mutations per locus.
    mutate_parent2();
    // Using params_no swap means that the following mutations will be passed from node zero
    // to node 2N:
    // 0, 0.51, 1, 2, 2.51, 3.0
    // Mutations at the following positions must be simplified out of the mutation table:
    // 1.51, 2.51
    // This, we have the following expectations:
    std::vector<double> muts_on_node_0
        = { 0., 0.51, 1., 2., 2.51, 3. }, //node_0 means post-simplification!
        muts_lost = { 1.51, 3.51 };
    BOOST_REQUIRE_EQUAL(tables.mutation_table.size(), 8);
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
    samples = { 0, 1 };
    fwdpp::ts::count_mutations(tables, pop.mutations, samples, pop.mcounts);
    decltype(pop.mcounts) mc(pop.mcounts.size(), 0);
    // NOTE a trick here:
    // The mutations transmitted from the parents are neutral.
    // I should probably fix that... ;)
    for (auto& locus : offspring)
        {
            for (auto m : pop.gametes[locus.first].smutations)
                {
                    mc[m]++;
                }
            for (auto m : pop.gametes[locus.first].mutations)
                {
                    mc[m]++;
                }
            for (auto m : pop.gametes[locus.second].smutations)
                {
                    mc[m]++;
                }
            for (auto m : pop.gametes[locus.second].mutations)
                {
                    mc[m]++;
                }
        }
    BOOST_REQUIRE_EQUAL(mc == pop.mcounts, true);
}

/* Tests based on a second, different set of w/in and b/w locus breakpoints
 *
 * The setup is as follows:
 *
 * Locus w/in b/w
 * 0     1    N/A
 * 1     2    1
 * 2     1    1
 * 3     0    0
 *
 * The expected breakpoints are as in expected_breakpoints in the fixture.
 */

BOOST_FIXTURE_TEST_CASE(test_breakpoints2, multilocus_fixture_deterministic)
{
    poptype::diploid_t offspring;

    BOOST_TEST_PASSPOINT();
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap2, offspring);
    BOOST_REQUIRE_EQUAL(data_to_record.first.swapped, 0);
    BOOST_REQUIRE_EQUAL(data_to_record.second.swapped, 0);
    auto expected = begin(expected_breakpoints2);
    auto brk = begin(data_to_record.first.breakpoints);

    BOOST_TEST_PASSPOINT();

    BOOST_REQUIRE_EQUAL(
        std::distance(expected, end(expected_breakpoints2)),
        std::distance(brk, end(data_to_record.first.breakpoints)));

    while (expected < end(expected_breakpoints2))
        {
            BOOST_REQUIRE_EQUAL(*expected, *brk);
            ++expected;
            ++brk;
        }
}

BOOST_FIXTURE_TEST_CASE(test_transmission2, multilocus_fixture_deterministic)
{
    mutate_parent2();
    poptype::diploid_t offspring;

    BOOST_TEST_PASSPOINT();
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap2, offspring);

    BOOST_TEST_PASSPOINT();
    BOOST_REQUIRE_NO_THROW(
        validate_mutations_positions_1_recparams2(offspring));

    //Detailed test of offspring gametes
    std::vector<double> expected_lost_muts = { 0.51, 1.51, 2. };
    std::vector<fwdpp::uint_t> keys_first, keys_second;
    for (auto& l : offspring)
        {
            for (auto m : pop.gametes[l.first].mutations)
                {
                    keys_first.push_back(m);
                }
            for (auto m : pop.gametes[l.second].mutations)
                {
                    keys_second.push_back(m);
                }
        }
    BOOST_REQUIRE_EQUAL(keys_second.size(), 0);
    for (auto p : expected_lost_muts)
        {
            auto itr = std::find_if(begin(keys_first), end(keys_first),
                                    [this, p](fwdpp::uint_t k) {
                                        return pop.mutations[k].pos == p;
                                    });
            BOOST_REQUIRE_EQUAL((itr == end(keys_first)), true);
            itr = std::find_if(begin(keys_second), end(keys_second),
                               [this, p](fwdpp::uint_t k) {
                                   return pop.mutations[k].pos == p;
                               });
            BOOST_REQUIRE_EQUAL((itr == end(keys_second)), true);
        }
}

BOOST_FIXTURE_TEST_CASE(test_edges2, multilocus_fixture_deterministic)
{
    poptype::diploid_t offspring;

    BOOST_TEST_PASSPOINT();
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap2, offspring);

    BOOST_TEST_PASSPOINT();

    fwdpp::ts::TS_NODE_INT next_index = tables.node_table.size();
    auto p1d = fwdpp::ts::get_parent_ids(0, 0, data_to_record.first.swapped);
    auto p2d = fwdpp::ts::get_parent_ids(0, 1, data_to_record.second.swapped);
    tables.add_offspring_data(next_index++, data_to_record.first.breakpoints,
                              data_to_record.first.mutation_keys, p1d, 0, 1);
    tables.add_offspring_data(next_index++, data_to_record.second.breakpoints,
                              data_to_record.second.mutation_keys, p2d, 0, 1);
    BOOST_REQUIRE_EQUAL(tables.edge_table.size(), 14);
    std::vector<std::pair<double, double>> expected_left_right
        = { { 0., 0.5 },  { 0.5, 1. }, { 1., 1.25 }, { 1.25, 1.75 },
            { 1.75, 2. }, { 2., 2.5 }, { 2.5, 4. } };
    std::size_t idx = 0, pidx = 0;
    fwdpp::ts::TS_NODE_INT parents[2] = { p1d.first, p1d.second };
    for (auto e : tables.edge_table)
        {
            BOOST_REQUIRE_EQUAL(e.left, expected_left_right[idx].first);
            BOOST_REQUIRE_EQUAL(e.right, expected_left_right[idx].second);
            BOOST_REQUIRE_EQUAL(parents[pidx], e.parent);
            ++idx;
            pidx = !pidx;
            if (idx == expected_left_right.size())
                {
                    idx = 0;
                    parents[0] = p2d.first;
                    parents[1] = p2d.second;
                    pidx = 0;
                }
        }
}

// Now, check with BOTH parents being mutated.
BOOST_FIXTURE_TEST_CASE(test_transmission_rec2_both_parents_mutated,
                        multilocus_fixture_deterministic)
{
    BOOST_REQUIRE_NO_THROW(mutate_both_parents());
    poptype::diploid_t offspring;

    BOOST_TEST_PASSPOINT();
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap2, offspring);

    BOOST_REQUIRE_EQUAL(
        data_to_record.first.breakpoints == expected_breakpoints2, true);
    BOOST_REQUIRE_EQUAL(
        data_to_record.second.breakpoints == expected_breakpoints2, true);

    std::vector<fwdpp::uint_t> g1keys, g2keys, g1keys_sel, g2keys_sel;
    for (auto& locus : offspring)
        {
            g1keys.insert(end(g1keys),
                          begin(pop.gametes[locus.first].mutations),
                          end(pop.gametes[locus.first].mutations));
            g2keys.insert(end(g2keys),
                          begin(pop.gametes[locus.second].mutations),
                          end(pop.gametes[locus.second].mutations));
            g1keys_sel.insert(end(g1keys_sel),
                              begin(pop.gametes[locus.first].smutations),
                              end(pop.gametes[locus.first].smutations));
            g2keys_sel.insert(end(g2keys_sel),
                              begin(pop.gametes[locus.second].smutations),
                              end(pop.gametes[locus.second].smutations));
        }
    // Validate gamete 1 of offspring
    for (auto p : expected_transmitted_mutations_mutate_both_parents_gamete_1)
        {
            auto itr = std::find_if(
                begin(g1keys), end(g1keys), [this, p](fwdpp::uint_t k) {
                    return std::fabs(pop.mutations[k].pos - p) <= 1e-6;
                });
            BOOST_REQUIRE_EQUAL(itr == end(g1keys), false);
        }
    // Validate gamete 2 of offspring
    for (auto p : expected_transmitted_mutations_mutate_both_parents_gamete_2)
        {
            auto itr = std::find_if(
                begin(g2keys), end(g2keys), [this, p](fwdpp::uint_t k) {
                    return std::fabs(pop.mutations[k].pos - p) <= 1e-6;
                });
            BOOST_REQUIRE_EQUAL(itr == end(g2keys), false);
        }

    // Finally, test that the RAW NUMBER of mutations inherited from time 0
    // == the expected number.  This is partly a test to ensure that
    // we've set the fixture up properly.
    int g1t0 = 0, g1t1 = 0, g2t0 = 0, g2t1 = 0;
    for (auto k : g1keys)
        {
            if (pop.mutations[k].g == 0)
                {
                    ++g1t0;
                }
        }
    for (auto k : g1keys_sel)
        {
            if (pop.mutations[k].g == 1)
                {
                    ++g1t1;
                }
        }
    for (auto k : g2keys)
        {
            if (pop.mutations[k].g == 0)
                {
                    ++g2t0;
                }
        }
    for (auto k : g2keys_sel)
        {
            if (pop.mutations[k].g == 1)
                {
                    ++g2t1;
                }
        }
    BOOST_REQUIRE_EQUAL(
        g1t0,
        expected_transmitted_mutations_mutate_both_parents_gamete_1.size());
    BOOST_REQUIRE_EQUAL(
        g2t0,
        expected_transmitted_mutations_mutate_both_parents_gamete_2.size());
    BOOST_REQUIRE_EQUAL(g1t1, 4);
    BOOST_REQUIRE_EQUAL(g2t1, 4);
}

BOOST_FIXTURE_TEST_CASE(test_ts_recording_rec2_both_parents_mutated,
                        multilocus_fixture_deterministic)
// Same data as previous test, but now we simplify.
{
    BOOST_REQUIRE_NO_THROW(mutate_both_parents());
    poptype::diploid_t offspring;

    BOOST_TEST_PASSPOINT();
    auto data_to_record = fwdpp::ts::generate_offspring(
        rng.get(), std::make_pair(0, 1), fwdpp::ts::selected_variants_only(),
        pop, params_no_swap2, offspring);
    fwdpp::ts::TS_NODE_INT next_index = tables.node_table.size();
    auto p1d = fwdpp::ts::get_parent_ids(0, 0, data_to_record.first.swapped);
    auto p2d = fwdpp::ts::get_parent_ids(0, 1, data_to_record.second.swapped);
    tables.add_offspring_data(next_index++, data_to_record.first.breakpoints,
                              data_to_record.first.mutation_keys, p1d, 0, 1);
    tables.add_offspring_data(next_index++, data_to_record.second.breakpoints,
                              data_to_record.second.mutation_keys, p2d, 0, 1);
    tables.sort_tables(pop.mutations);

    // All expected variants must be in the mutation table PRIOR TO SIMPLIFICATION
    for (auto p : expected_transmitted_mutations_mutate_both_parents_gamete_1)
        {
            bool found = false;
            for (auto m : tables.mutation_table)
                {
                    if (std::fabs(pop.mutations[m.key].pos - p) <= 1e-6)
                        {
                            found = true;
                            break;
                        }
                }
            BOOST_REQUIRE_EQUAL(found, true);
        }

    for (auto p : expected_transmitted_mutations_mutate_both_parents_gamete_2)
        {
            bool found = false;
            for (auto m : tables.mutation_table)
                {
                    if (std::fabs(pop.mutations[m.key].pos - p) <= 1e-6)
                        {
                            found = true;
                            break;
                        }
                }
            BOOST_REQUIRE_EQUAL(found, true);
        }
    fwdpp::ts::table_simplifier simplifier(nloci);
    std::vector<fwdpp::ts::TS_NODE_INT> samples(
        { next_index - 2, next_index - 1 });
    for(auto n : samples){BOOST_REQUIRE_EQUAL(tables.node_table[n].time,1.);}
    auto rv = simplifier.simplify(tables, samples, pop.mutations);

    // All expected variants must be in the mutation table!
    for (auto p : expected_transmitted_mutations_mutate_both_parents_gamete_1)
        {
            bool found = false;
            for (auto m : tables.mutation_table)
                {
                    if (std::fabs(pop.mutations[m.key].pos - p) <= 1e-6)
                        {
                            found = true;
                            break;
                        }
                }
            BOOST_CHECK_EQUAL(found, true);
        }

    for (auto p : expected_transmitted_mutations_mutate_both_parents_gamete_2)
        {
            bool found = false;
            for (auto m : tables.mutation_table)
                {
                    if (std::fabs(pop.mutations[m.key].pos - p) <= 1e-6)
                        {
                            found = true;
                            break;
                        }
                }
            BOOST_TEST_CHECKPOINT("Checking for mutation in gamete 2 after "
                                  "simplification at position "
                                  << p);
            BOOST_CHECK_EQUAL(found, true);
        }
    // Count the number of mutations surviving simplification
    // from generation 0 to 1 (parents to offspring):

    int muts_sample_0 = 0, muts_sample_1 = 0;
    int new_muts_0 = 0, new_muts_1 = 0;
    for (auto& m : tables.mutation_table)
        {
            if (m.node == 0 && pop.mutations[m.key].g == 0)
                {
                    ++muts_sample_0;
                }
            if (m.node == 0 && pop.mutations[m.key].g == 1)
                {
                    ++new_muts_0;
                }
            if (m.node == 1 && pop.mutations[m.key].g == 0)
                {
                    ++muts_sample_1;
                }
            if (m.node == 1 && pop.mutations[m.key].g == 1)
                {
                    ++new_muts_1;
                }
        }
    BOOST_REQUIRE_EQUAL(new_muts_0,4);
    BOOST_REQUIRE_EQUAL(new_muts_1,4);
    BOOST_REQUIRE_EQUAL(
        muts_sample_0,
        expected_transmitted_mutations_mutate_both_parents_gamete_1.size());
    BOOST_REQUIRE_EQUAL(
        muts_sample_1,
        expected_transmitted_mutations_mutate_both_parents_gamete_2.size());
}
