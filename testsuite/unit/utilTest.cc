/*!
  \file utilTest.cc
  \ingroup unit
  \brief Testing that functions in fwdpp/util.hpp work
*/

#include <config.h>
#include <boost/test/unit_test.hpp>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/util.hpp>

/*
  A note on test names:
  All recycling is applied to extinct mutations.
  What varies from test to test is the treatment of fixations.
*/

BOOST_AUTO_TEST_SUITE(test_util)

BOOST_FIXTURE_TEST_CASE(only_recycle_extinct_mutations,
                        standard_empty_single_deme_fixture)
{
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    mcounts.emplace_back(1);
    mutations.emplace_back(0.2, 0);
    mut_lookup.emplace(0.2, 1);
    mcounts.emplace_back(0);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }
    fwdpp::update_mutations(mutations, mut_lookup, mcounts);
    BOOST_REQUIRE_EQUAL(mcounts[0], 1);
    BOOST_REQUIRE_EQUAL(mcounts[1], 0);
    BOOST_REQUIRE_EQUAL(mutations.size(), 2);
    BOOST_REQUIRE_EQUAL(fixations.size(), 0);
    BOOST_REQUIRE_EQUAL(fixation_times.size(), 0);
}

BOOST_FIXTURE_TEST_CASE(test_recycle_all_fixations,
                        standard_empty_single_deme_fixture)
{
    fwdpp::uint_t N = 1000;
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }
    mcounts.emplace_back(2 * N);
    fwdpp::update_mutations(mutations, mut_lookup, mcounts, 2 * N);
    BOOST_REQUIRE_EQUAL(mcounts[0], 0);
    BOOST_REQUIRE_EQUAL(mutations.size(), 1);
    BOOST_REQUIRE_EQUAL(fixations.size(), 0);
    BOOST_REQUIRE_EQUAL(fixation_times.size(), 0);
}

BOOST_FIXTURE_TEST_CASE(
    test_recycle_all_fixations_and_move_to_fixations_vector,
    standard_empty_single_deme_fixture)
{
    fwdpp::uint_t N = 1000;
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }
    mcounts.emplace_back(2 * N);
    // use generation = 2
    fwdpp::update_mutations(mutations, fixations, fixation_times, mut_lookup,
                            mcounts, 2, 2 * N);
    BOOST_REQUIRE_EQUAL(mcounts[0], 0);
    BOOST_REQUIRE_EQUAL(mutations.size(), 1);
    BOOST_REQUIRE_EQUAL(fixations.size(), 1);
    BOOST_REQUIRE_EQUAL(fixation_times.size(), 1);
    BOOST_REQUIRE_EQUAL(fixation_times[0], 2);
}

BOOST_FIXTURE_TEST_CASE(only_recycle_neutral_fixations,
                        standard_empty_single_deme_fixture)
{
    fwdpp::uint_t N = 1000;
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    mcounts.emplace_back(2 * N);
    // This is NOT a neutral mutation
    mutations.emplace_back(0.2, -1.0);
    mut_lookup.emplace(0.2, 1);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }
    mcounts.emplace_back(2 * N);
    // use generation = 2
    fwdpp::update_mutations_n(mutations, fixations, fixation_times, mut_lookup,
                              mcounts, 2, 2 * N);
    BOOST_REQUIRE_EQUAL(mcounts[0], 0);
    BOOST_REQUIRE_EQUAL(mcounts[1], 2 * N);
    BOOST_REQUIRE_EQUAL(mutations.size(), 2);
    BOOST_REQUIRE_EQUAL(fixations.size(), 2);
    BOOST_REQUIRE_EQUAL(fixation_times.size(), 2);
    BOOST_REQUIRE_EQUAL(fixation_times[0], 2);
    // Do it again, but
    // use generation = 3
    fwdpp::update_mutations_n(mutations, fixations, fixation_times, mut_lookup,
                              mcounts, 3, 2 * N);
    // This will not change the data b/c
    // the neutral variant is already flagged to recycle
    // and the selected one is already entered into fixations.
    BOOST_REQUIRE_EQUAL(mcounts[0], 0);
    BOOST_REQUIRE_EQUAL(mcounts[1], 2 * N);
    BOOST_REQUIRE_EQUAL(mutations.size(), 2);
    BOOST_REQUIRE_EQUAL(fixations.size(), 2);
    BOOST_REQUIRE_EQUAL(fixation_times.size(), 2);
    BOOST_REQUIRE_EQUAL(fixation_times[0], 2);
    BOOST_REQUIRE_EQUAL(fixation_times[1], 2);
}

// The following unit tests are checks for GitHub issue #130
// This issue involves the situation where an extinct mutation
// has the same position as an extant mutation.  Through fwdpp
// 0.6.0, this cause the extant mutation's position to
// be removed from the lookup table, which is a book-keeping bug.
// However, I ran a ton of sims and never found a difference due
// to the bug vs the fixed version (fwdpp 0.6.1).  The reason
// is that, for the bug to have an effect, a series of rare events
// must happen multiple times.

BOOST_FIXTURE_TEST_CASE(issue_130_diploid_population_test_update_mutations,
                        standard_empty_single_deme_fixture)
{
    fwdpp::uint_t N = 1000;
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    mcounts.emplace_back(2 * N);
    // This is NOT a neutral mutation
    mutations.emplace_back(0.2, -1.0);
    mut_lookup.emplace(0.2, 1);
    mcounts.push_back(1); //a singleton

    // This is an exctinct variant at 0.2
    mutations.emplace_back(0.2, 0.0);
    mcounts.push_back(0);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }

    fwdpp::update_mutations(mutations, fixations, fixation_times, mut_lookup,
                            mcounts, 0, 2 * N);
    for (std::size_t i = 0; i < mcounts.size(); ++i)
        {
            if (mcounts[i] > 0)
                {
                    BOOST_CHECK_EQUAL(mut_lookup.find(mutations[i].pos)
                                          != mut_lookup.end(),
                                      true);
                }
        }
}

BOOST_FIXTURE_TEST_CASE(issue_130_diploid_population_test_update_mutations_extinct_only,
                        standard_empty_single_deme_fixture)
{
    fwdpp::uint_t N = 1000;
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    mcounts.emplace_back(2 * N);
    // This is NOT a neutral mutation
    mutations.emplace_back(0.2, -1.0);
    mut_lookup.emplace(0.2, 1);
    mcounts.push_back(1); //a singleton

    // This is an exctinct variant at 0.2
    mutations.emplace_back(0.2, 0.0);
    mcounts.push_back(0);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }

    fwdpp::update_mutations(mutations, mut_lookup, mcounts);
    for (std::size_t i = 0; i < mcounts.size(); ++i)
        {
            if (mcounts[i] > 0)
                {
                    BOOST_CHECK_EQUAL(mut_lookup.find(mutations[i].pos)
                                          != mut_lookup.end(),
                                      true);
                }
        }
}
BOOST_FIXTURE_TEST_CASE(
    issue_130_diploid_population_test_update_mutations_do_not_record_fixations,
    standard_empty_single_deme_fixture)
{
    fwdpp::uint_t N = 1000;
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    mcounts.emplace_back(2 * N);
    // This is NOT a neutral mutation
    mutations.emplace_back(0.2, -1.0);
    mut_lookup.emplace(0.2, 1);
    mcounts.push_back(1); //a singleton

    // This is an exctinct variant at 0.2
    mutations.emplace_back(0.2, 0.0);
    mcounts.push_back(0);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }

    fwdpp::update_mutations(mutations, mut_lookup, mcounts, 2 * N);
    for (std::size_t i = 0; i < mcounts.size(); ++i)
        {
            if (mcounts[i] > 0)
                {
                    BOOST_CHECK_EQUAL(mut_lookup.find(mutations[i].pos)
                                          != mut_lookup.end(),
                                      true);
                }
        }
}

BOOST_FIXTURE_TEST_CASE(issue_130_diploid_population_test_update_mutations_n,
                        standard_empty_single_deme_fixture)
{
    fwdpp::uint_t N = 1000;
    // This is a neutral mutation
    mutations.emplace_back(0.1, 0);
    mut_lookup.emplace(0.1, 0);
    mcounts.emplace_back(2 * N);
    // This is NOT a neutral mutation
    mutations.emplace_back(0.2, -1.0);
    mut_lookup.emplace(0.2, 1);
    mcounts.push_back(1); //a singleton

    // This is an exctinct variant at 0.2
    mutations.emplace_back(0.2, 0.0);
    mcounts.push_back(0);
    for (auto& i : mut_lookup)
        {
            BOOST_REQUIRE_EQUAL(i.first, mutations[i.second].pos);
        }

    fwdpp::update_mutations_n(mutations, fixations, fixation_times, mut_lookup,
                              mcounts, 0, 2 * N);
    for (std::size_t i = 0; i < mcounts.size(); ++i)
        {
            if (mcounts[i] > 0)
                {
                    BOOST_CHECK_EQUAL(mut_lookup.find(mutations[i].pos)
                                          != mut_lookup.end(),
                                      true);
                }
        }
}

BOOST_AUTO_TEST_SUITE_END()
