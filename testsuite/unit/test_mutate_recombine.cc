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
    std::vector<double> breakpoints(1, 1.5);
    auto q = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);
    BOOST_REQUIRE_EQUAL(q.empty(), true);
    

    auto new_gamete_index = fwdpp::mutate_recombine(
        new_mutations, breakpoints, 0, 1, pop.gametes, pop.mutations, q,
        pop.neutral, pop.selected);

    // The test is that the new gamete does NOT contain a mutation
    // at pos 1.5
    auto itr = std::find(pop.gametes[new_gamete_index].mutations.begin(),
                         pop.gametes[new_gamete_index].mutations.end(),
                         pop.mutations.size() - 1);
    auto is_end = (itr == pop.gametes[new_gamete_index].mutations.end());
    BOOST_REQUIRE_EQUAL(is_end, true);

    // But, if we switch the parent IDs,
    // then this mutation WILL appear in the offspring gamete
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
}

BOOST_AUTO_TEST_SUITE_END()
