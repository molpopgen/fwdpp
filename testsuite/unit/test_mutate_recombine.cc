#include <iostream>
#include <cassert>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>

BOOST_AUTO_TEST_SUITE(test_mutate_recombine)

BOOST_FIXTURE_TEST_CASE(test_boundary, slocuspop_objects)
{
    poptype pop(diploids, gametes, mutations);

    // We are going to manually add a neutral mutation at position 1.5
    pop.mutations.emplace_back(1.5, 0.0, 1.0, 1);
    pop.mcounts.push_back(1); // mock the mutation count
    // We add the mutation to the second gamete of the pop
    fwdpp::fwdpp_internal::insert_new_mutation(
        pop.gametes[1].mutations.begin(), pop.gametes[1].mutations.end(),
        pop.mutations.size() - 1, pop.mutations, pop.gametes[1].mutations);

    BOOST_REQUIRE_EQUAL(pop.mutations.back().pos, 1.5);
    std::vector<fwdpp::uint_t> new_mutations;
    std::vector<double> breakpoints{ 1.5, std::numeric_limits<double>::max() };
    auto q = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);
    BOOST_REQUIRE_EQUAL(q.empty(), true);

    auto new_gamete_index = fwdpp::mutate_recombine(
        new_mutations, breakpoints, 0, 1, pop.gametes, pop.mutations, q,
        pop.neutral, pop.selected);

    // The test is that the new gamete DOES contain a mutation
    // at pos 1.5
    auto itr = std::find(pop.gametes[new_gamete_index].mutations.begin(),
                         pop.gametes[new_gamete_index].mutations.end(),
                         pop.mutations.size() - 1);
    auto is_end = (itr == pop.gametes[new_gamete_index].mutations.end());
    BOOST_REQUIRE_EQUAL(is_end, false);

    // But, if we switch the parent IDs,
    // then this mutation WILL appear in the offspring gamete
    new_gamete_index = fwdpp::mutate_recombine(new_mutations, breakpoints, 1,
                                               0, pop.gametes, pop.mutations,
                                               q, pop.neutral, pop.selected);

    // The test is that the new gamete DOES NOT contain a mutation
    // at pos 1.5
    itr = std::find(pop.gametes[new_gamete_index].mutations.begin(),
                    pop.gametes[new_gamete_index].mutations.end(),
                    pop.mutations.size() - 1);
    is_end = (itr == pop.gametes[new_gamete_index].mutations.end());
    BOOST_REQUIRE_EQUAL(is_end, true);
}

BOOST_FIXTURE_TEST_CASE(test_boundary_with_recurrent_mutation,
                        slocuspop_objects)
    /// This test is similar to the above, but we use the xtra field of
    /// fwdpp::popgenmut for extra resolution about how mutation vs
    /// inherited mutations are treated vis-a-vis crossover breakpoints
{
    poptype pop(diploids, gametes, mutations);

    // We are going to manually add a neutral mutation at position 1.5
    pop.mutations.emplace_back(1.5, 0.0, 1.0, 1, 7);
    pop.mcounts.push_back(1); // mock the mutation count
    // We add the mutation to the second gamete of the pop
    fwdpp::fwdpp_internal::insert_new_mutation(
        pop.gametes[1].mutations.begin(), pop.gametes[1].mutations.end(),
        pop.mutations.size() - 1, pop.mutations, pop.gametes[1].mutations);

    BOOST_REQUIRE_EQUAL(pop.mutations.back().pos, 1.5);
    // We will now add another mutation at 1.5, which will be
    // a new mutation in the offspring
    pop.mutations.emplace_back(1.5, 0.0, 1.0, 1, 9);
    BOOST_REQUIRE_EQUAL(pop.mutations.back().pos, 1.5);
    // The new mutation at 1.5 is at the end
    std::vector<fwdpp::uint_t> new_mutations(1, pop.mutations.size() - 1);
    std::vector<double> breakpoints{ 1.5, std::numeric_limits<double>::max() };
    auto q = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);
    BOOST_REQUIRE_EQUAL(q.empty(), true);

    auto new_gamete_index = fwdpp::mutate_recombine(
        new_mutations, breakpoints, 0, 1, pop.gametes, pop.mutations, q,
        pop.neutral, pop.selected);

    // The test is that the new gamete DOES contain a mutation
    // at pos 1.5
    auto itr = std::find(pop.gametes[new_gamete_index].mutations.begin(),
                         pop.gametes[new_gamete_index].mutations.end(),
                         pop.mutations.size() - 1);
    auto is_end = (itr == pop.gametes[new_gamete_index].mutations.end());
    BOOST_REQUIRE_EQUAL(is_end, false);

    //Further, there must be two mutations at position 1.5, because one
    //is inherited from the second parental gamete and the other is
    //a new mutation.
    auto count_new_mutation
        = std::count_if(pop.gametes[new_gamete_index].mutations.begin(),
                        pop.gametes[new_gamete_index].mutations.end(),
                        [&pop](const fwdpp::uint_t key) {
                            return pop.mutations[key].pos == 1.5;
                        });
    BOOST_REQUIRE_EQUAL(count_new_mutation, 2);

    // But, if we switch the parent IDs,
    // then this mutation WILL NOT appear in the offspring gamete
    new_gamete_index = fwdpp::mutate_recombine(new_mutations, breakpoints, 1,
                                               0, pop.gametes, pop.mutations,
                                               q, pop.neutral, pop.selected);

    // The test is that the new gamete DOES contain a mutation
    // at pos 1.5
    itr = std::find(pop.gametes[new_gamete_index].mutations.begin(),
                    pop.gametes[new_gamete_index].mutations.end(),
                    pop.mutations.size() - 1);
    is_end = (itr == pop.gametes[new_gamete_index].mutations.end());
    BOOST_REQUIRE_EQUAL(is_end, false);
    // However, it must only contain one copy, and it must be the 
    // new mutation, whose xtra field is set to 9
    BOOST_REQUIRE_EQUAL(static_cast<int>(pop.mutations[*itr].xtra), 9);
    count_new_mutation
        = std::count_if(pop.gametes[new_gamete_index].mutations.begin(),
                        pop.gametes[new_gamete_index].mutations.end(),
                        [&pop](const fwdpp::uint_t key) {
                            return pop.mutations[key].pos == 1.5;
                        });
    BOOST_REQUIRE_EQUAL(count_new_mutation, 1);
}

BOOST_AUTO_TEST_SUITE_END()
