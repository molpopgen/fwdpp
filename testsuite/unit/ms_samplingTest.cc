/*!
  \file ms_sampling.cc
  \ingroup unit
  \brief Testing KTfwd::ms_sample
*/

#include <iostream>
#include <config.h>
#include <boost/test/unit_test.hpp>
#include <fwdpp/diploid.hh>
#include "../fixtures/fwdpp_fixtures.hpp"
#include <fwdpp/sampling_functions.hpp>

BOOST_FIXTURE_TEST_SUITE(ms_samplingTest, standard_empty_single_deme_fixture)

void
add_mutations(mcont_t &mutations, gcont_t &gametes, dipvector_t &diploids)
{
    // auto first_gam = gametes.begin();
    for (unsigned i = 0; i < 10; ++i)
        {
            mutations.emplace_back(double(i), 0., 0.);
            gcont_t::value_type ng(1);
            ng.mutations.push_back(mutations.size() - 1);

            mutations.emplace_back(double(i) + 0.5, 0., 0.);
            gcont_t::value_type ng2(1);
            ng2.mutations.push_back(mutations.size() - 1);

            gametes.emplace_back(std::move(ng));
            diploids[i].first = gametes.size() - 1;
            gametes.emplace_back(std::move(ng2));
            diploids[i].second = gametes.size() - 1;
        }
    gametes[0].n = 0;
    // gametes.erase(first_gam);
}

BOOST_AUTO_TEST_CASE(sample_properties_test)
{
    diploids.resize(10, std::make_pair(0, 0));
    gametes.emplace_back(20);
    add_mutations(mutations, gametes, diploids);
    auto sample = KTfwd::fwdpp_internal::ms_sample_separate_single_deme(
        mutations, gametes, diploids,
        std::vector<unsigned>({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }), 20, true);
    BOOST_REQUIRE_EQUAL(sample.second.empty(), true); // no selected mutations
    BOOST_REQUIRE_EQUAL(
        std::is_sorted(sample.first.begin(), sample.first.end(),
                       [](const std::pair<double, std::string> &a,
                          const std::pair<double, std::string> &b) {
                           return a.first < b.first;
                       }),
        true);
    for (unsigned i = 0; i < sample.first.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(sample.first[i].second.size(), 20);
        }
    BOOST_REQUIRE_EQUAL(sample.first.size(), 20); // 20 neutral mutations
}

BOOST_AUTO_TEST_CASE(odd)
{
    diploids.resize(10, std::make_pair(0, 0));
    gametes.emplace_back(20);
    add_mutations(mutations, gametes, diploids);
    auto sample = KTfwd::fwdpp_internal::ms_sample_separate_single_deme(
        mutations, gametes, diploids,
        std::vector<unsigned>({ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }), 19, true);
    BOOST_REQUIRE_EQUAL(sample.second.empty(), true); // no selected mutations
    BOOST_REQUIRE_EQUAL(
        std::is_sorted(sample.first.begin(), sample.first.end(),
                       [](const std::pair<double, std::string> &a,
                          const std::pair<double, std::string> &b) {
                           return a.first < b.first;
                       }),
        true);
    for (unsigned i = 0; i < sample.first.size(); ++i)
        {
            BOOST_REQUIRE_EQUAL(sample.first[i].second.size(), 19);
        }
    // Aha--need to have a "remove invariant"
    BOOST_REQUIRE_EQUAL(sample.first.size(), 19); // 20 neutral mutations
}

BOOST_AUTO_TEST_SUITE_END()
