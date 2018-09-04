// Test the construction of populations
// from user input data

#include <cassert>
#include <numeric>
#include <boost/test/unit_test.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>
#include <testsuite/util/quick_evolve_sugar.hpp>
#include <gsl/gsl_matrix_char.h>

BOOST_FIXTURE_TEST_SUITE(test_sampling, slocuspop_objects)

BOOST_AUTO_TEST_CASE(test_taking_sample)
{
    diploids.emplace_back(0, 1);
    gametes[0].n++;
    gametes[1].n++;
    poptype p(diploids, gametes, mutations);
    auto m = fwdpp::sample_individuals(p, std::vector<std::size_t>{ 0, 1 },
                                       true, true, true);
    for (auto k : m.neutral_keys)
        {
            BOOST_REQUIRE_EQUAL(p.mutations[k].neutral, true);
        }
    for (auto k : m.selected_keys)
        {
            BOOST_REQUIRE_EQUAL(p.mutations[k].neutral, false);
        }
    BOOST_REQUIRE_EQUAL(m.neutral_keys.size(), 1);
    BOOST_REQUIRE_EQUAL(m.selected_keys.size(), 2);
    BOOST_REQUIRE_EQUAL(m.neutral.data.size(), 4);
    BOOST_REQUIRE_EQUAL(m.selected.data.size(), 8);
    BOOST_REQUIRE_EQUAL(
        std::accumulate(m.neutral.data.begin(), m.neutral.data.end(), 0), 2);
    BOOST_REQUIRE_EQUAL(
        std::accumulate(m.selected.data.begin(), m.selected.data.end(), 0), 4);
    for (std::size_t i = 0; i < m.neutral.data.size(); ++i)
        {
            if (i % 2 == 0.0)
                {
                    BOOST_REQUIRE_EQUAL(m.neutral.data[i], 1);
                }
            else
                {
                    BOOST_REQUIRE_EQUAL(m.neutral.data[i], 0);
                }
        }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_sampling_mlocus, mlocuspop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(test_windowed_sampling)
{
    simulate_mlocuspop(pop, rng, mutmodels, recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), mu,
                       rbw, generation, 1000);
    BOOST_CHECK_EQUAL(generation, 1000);
    std::vector<std::size_t> individuals;
    for (unsigned i = 0; i < 50; ++i)
        {
            auto x = static_cast<std::size_t>(
                gsl_ran_flat(rng.get(), 0, pop.diploids.size()));
            while (std::find(individuals.begin(), individuals.end(), x)
                   != individuals.end())
                {
                    x = static_cast<std::size_t>(
                        gsl_ran_flat(rng.get(), 0, pop.diploids.size()));
                }
            individuals.push_back(x);
        }
    auto vdm = fwdpp::sample_individuals_by_window(
        pop, individuals, pop.locus_boundaries, true, false, true);
    BOOST_REQUIRE_EQUAL(vdm.size(), pop.locus_boundaries.size());
    // To test contents, we'll pull the mutation keys that were used
    auto keys = fwdpp::fwdpp_internal::generate_filter_sort_keys(
        pop, individuals, true, false, true);
    unsigned totmuts = 0;
    for (auto& i : vdm)
        {
            totmuts += i.neutral.positions.size();
        }
    BOOST_REQUIRE_EQUAL(totmuts, keys.first.size());

    // The more interesting test is to ensure that each locus is correct,
    // which we can brute-force:

    for (std::size_t i = 0; i < vdm.size(); ++i)
        {
            auto m = gsl_matrix_char_const_view_array(
                reinterpret_cast<const char*>(vdm[i].neutral.data.data()),
                vdm[i].neutral.positions.size(), 2 * individuals.size());
            for (std::size_t row = 0; row < m.matrix.size1; ++row)
                {
                    int rsum = 0;
                    for (std::size_t hap = 0; hap < m.matrix.size2; ++hap)
                        {
                            std::int8_t state(static_cast<std::int8_t>(
                                gsl_matrix_char_get(&m.matrix, row, hap)));
                            rsum += state;
                            std::size_t ind_index = hap / 2;
                            BOOST_REQUIRE_EQUAL(ind_index < individuals.size(),
                                                true);
                            auto ind = individuals[ind_index];
                            int gam = (hap % 2 == 0) ? 0 : 1;
                            auto pos = vdm[i].neutral.positions[row];
                            const auto& gamref
                                = (gam == 0)
                                      ? pop.gametes[pop.diploids[ind][i].first]
                                      : pop.gametes[pop.diploids[ind][i]
                                                        .second];
                            auto itr = std::find_if(
                                gamref.mutations.begin(),
                                gamref.mutations.end(),
                                [pos, this](const fwdpp::uint_t k) {
                                    return pop.mutations[k].pos == pos;
                                });
                            int found
                                = (itr == gamref.mutations.end()) ? 0 : 1;
                            if (state == 0)
                                {
                                    BOOST_REQUIRE_EQUAL(found, 0);
                                }
                            else
                                {
                                    BOOST_REQUIRE_EQUAL(found, 1);
                                }
                        }
                    //Get the row sum directly from the data_matrix
                    auto msum
                        = std::accumulate(vdm[i].neutral.data.begin()
                                              + row * 2 * individuals.size(),
                                          vdm[i].neutral.data.begin()
                                              + row * 2 * individuals.size()
                                              + 2 * individuals.size(),
                                          0);
                    BOOST_REQUIRE_EQUAL(rsum, msum);
                }
        }
}

BOOST_AUTO_TEST_CASE(test_sampling_mlocus_empty_locus)
{
    mu[1] = 0.0; //make one of the loci get zero mutations
    simulate_mlocuspop(pop, rng, mutmodels, recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(), mu,
                       rbw, generation, 1000);
    BOOST_CHECK_EQUAL(generation, 1000);
    std::vector<std::size_t> individuals;
    for (unsigned i = 0; i < 50; ++i)
        {
            auto x = static_cast<std::size_t>(
                gsl_ran_flat(rng.get(), 0, pop.diploids.size()));
            while (std::find(individuals.begin(), individuals.end(), x)
                   != individuals.end())
                {
                    x = static_cast<std::size_t>(
                        gsl_ran_flat(rng.get(), 0, pop.diploids.size()));
                }
            individuals.push_back(x);
        }
    auto vdm = fwdpp::sample_individuals_by_window(
        pop, individuals, pop.locus_boundaries, true, false, true);
    BOOST_REQUIRE_EQUAL(vdm.size(), pop.locus_boundaries.size());
    // To test contents, we'll pull the mutation keys that were used
    auto keys = fwdpp::fwdpp_internal::generate_filter_sort_keys(
        pop, individuals, true, false, true);
    unsigned totmuts = 0;
    for (auto& i : vdm)
        {
            totmuts += i.neutral.positions.size();
        }
    BOOST_REQUIRE_EQUAL(totmuts, keys.first.size());

    // The more interesting test is to ensure that each locus is correct,
    // which we can brute-force:
    BOOST_REQUIRE_EQUAL(vdm[1].neutral.data.empty(),true);
    BOOST_REQUIRE_EQUAL(vdm[1].neutral.positions.empty(),true);
    for (std::size_t i = 0; i < vdm.size(); ++i)
        {
            for (auto p : vdm[i].neutral.positions)
                {
                    BOOST_CHECK_EQUAL(p >= pop.locus_boundaries[i].first,
                                      true);
                    BOOST_CHECK_EQUAL(p < pop.locus_boundaries[i].second,
                                      true);
                }
            if (!vdm[i].neutral.positions.empty())
                {
                    BOOST_REQUIRE_EQUAL(vdm[i].neutral.data.empty(), false);
                    auto m = gsl_matrix_char_const_view_array(
                        reinterpret_cast<const char*>(
                            vdm[i].neutral.data.data()),
                        vdm[i].neutral.positions.size(),
                        2 * individuals.size());
                    for (std::size_t row = 0; row < m.matrix.size1; ++row)
                        {
                            int rsum = 0;
                            for (std::size_t hap = 0; hap < m.matrix.size2;
                                 ++hap)
                                {
                                    std::int8_t state(static_cast<std::int8_t>(
                                        gsl_matrix_char_get(&m.matrix, row,
                                                            hap)));
                                    rsum += state;
                                    std::size_t ind_index = hap / 2;
                                    BOOST_REQUIRE_EQUAL(
                                        ind_index < individuals.size(), true);
                                    auto ind = individuals[ind_index];
                                    int gam = (hap % 2 == 0) ? 0 : 1;
                                    auto pos = vdm[i].neutral.positions[row];
                                    const auto& gamref
                                        = (gam == 0)
                                              ? pop.gametes
                                                    [pop.diploids[ind][i]
                                                         .first]
                                              : pop.gametes
                                                    [pop.diploids[ind][i]
                                                         .second];
                                    auto itr = std::find_if(
                                        gamref.mutations.begin(),
                                        gamref.mutations.end(),
                                        [pos, this](const fwdpp::uint_t k) {
                                            return pop.mutations[k].pos == pos;
                                        });
                                    int found = (itr == gamref.mutations.end())
                                                    ? 0
                                                    : 1;
                                    if (state == 0)
                                        {
                                            BOOST_REQUIRE_EQUAL(found, 0);
                                        }
                                    else
                                        {
                                            BOOST_REQUIRE_EQUAL(found, 1);
                                        }
                                }
                            //Get the row sum directly from the data_matrix
                            auto msum = std::accumulate(
                                vdm[i].neutral.data.begin()
                                    + row * 2 * individuals.size(),
                                vdm[i].neutral.data.begin()
                                    + row * 2 * individuals.size()
                                    + 2 * individuals.size(),
                                0);
                            BOOST_REQUIRE_EQUAL(rsum, msum);
                        }
                }
        }
}
BOOST_AUTO_TEST_SUITE_END()
