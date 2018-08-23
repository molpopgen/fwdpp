/* \file sugar_matrixTest.cc
 * \brief Integration and unit tests of data matrix generation
 * \ingroup unit
 */

#include <config.h>
#include <boost/test/unit_test.hpp>
#include "../fixtures/sugar_fixtures.hpp"
#include "../util/quick_evolve_sugar.hpp"
#include <fwdpp/data_matrix.hpp>
#include <fwdpp/sugar/sampling.hpp>
#include <gsl/gsl_matrix_char.h>

// This is an involved integration test of
// haplotype and genotype matrices, and it takes
// some time to run.
// We simulate a population of N=1000 diploids to equilibrium.
// We then construct a haplotype and genotype matrix from a large
// sample of diploids.  We tests that row and sum columns from the two
// matrix types agree.  Once those agreement tests pass, we compare the
// haplotype matrix to the row/sum columns of a sample based on
// fwdpp::sample_separate (fwdpp/sugar/sampling.hpp) for the same diploids.
// These last checks ensure that two independent pieces of code, written
// at different times, give the same results.
//
// Then, to really make sure, we do the above tests 1,000 times, evolving
// our population 100 generations in between each test.
BOOST_AUTO_TEST_CASE(slocuspop_hapmatrix_exhaustive)
{
    using spoptype = slocuspop_popgenmut_fixture::poptype;
    spoptype pop(1000);
    fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2> rng(0u);
    auto generation = simulate_slocuspop(pop, rng, 0, 10000);
    std::vector<std::size_t> indlist;
    // Sample a LOT of individuals
    for (std::size_t i = 100; i < 750; i += 5)
        indlist.push_back(i);
    for (unsigned ntests = 0; ntests < 1000; ++ntests)
        {
            auto keys = mutation_keys(pop, indlist, true, true);
            // keys come out unsorted, so we have to sort for this test:
            std::sort(keys.first.begin(), keys.first.end(),
                      [&pop](const std::pair<std::size_t, fwdpp::uint_t> &a,
                             const std::pair<std::size_t, fwdpp::uint_t> &b) {
                          return pop.mutations[a.first].pos
                                 < pop.mutations[b.first].pos;
                      });
            std::sort(keys.second.begin(), keys.second.end(),
                      [&pop](const std::pair<std::size_t, fwdpp::uint_t> &a,
                             const std::pair<std::size_t, fwdpp::uint_t> &b) {
                          return pop.mutations[a.first].pos
                                 < pop.mutations[b.first].pos;
                      });
            // In order to be comparable to a sample taken from the population,
            // we must remove sites fixed in these inidivudals
            keys.first.erase(
                std::remove_if(
                    keys.first.begin(), keys.first.end(),
                    [&indlist](
                        const std::pair<std::size_t, fwdpp::uint_t> &a) {
                        return a.second == 2 * indlist.size();
                    }),
                keys.first.end());
            keys.second.erase(
                std::remove_if(
                    keys.second.begin(), keys.second.end(),
                    [&indlist](
                        const std::pair<std::size_t, fwdpp::uint_t> &a) {
                        return a.second == 2 * indlist.size();
                    }),
                keys.second.end());
            // Create data matrices and do basic sanity checks
            auto m = haplotype_matrix(pop, indlist, keys.first, keys.second);
            BOOST_REQUIRE_EQUAL(m.ncol, 2 * indlist.size());
            BOOST_REQUIRE_EQUAL(m.neutral.size(), m.ncol * keys.first.size());
            BOOST_REQUIRE_EQUAL(m.selected.size(),
                                m.ncol * keys.second.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(), m.neutral_keys.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(),
                                m.selected_positions.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(), m.selected_keys.size());

            // Note that we can convert the data to vectors of other types
            // as needed
            std::vector<double> neutral_as_double(m.neutral.begin(),
                                                  m.neutral.end());
            for (std::size_t i = 0; i < m.neutral.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(m.neutral[i], neutral_as_double[i]);
                }

            auto gm = genotype_matrix(pop, indlist, keys.first, keys.second);
            BOOST_REQUIRE_EQUAL(gm.ncol, indlist.size());
            BOOST_REQUIRE_EQUAL(gm.neutral.size(),
                                gm.ncol * keys.first.size());
            BOOST_REQUIRE_EQUAL(gm.selected.size(),
                                gm.ncol * keys.second.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(),
                                gm.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(keys.first.size(), gm.neutral_keys.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(),
                                gm.selected_positions.size());
            BOOST_REQUIRE_EQUAL(keys.second.size(),
                                gm.selected_keys.size());

            // Now, compare to an independent calculation of the genotypes
            // for same individuals
            auto s = fwdpp::sample_separate(pop, indlist, true);
            // Check same total # of mutations
            BOOST_REQUIRE_EQUAL(s.first.size(), m.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(s.second.size(), m.selected_positions.size());
            BOOST_REQUIRE_EQUAL(s.first.size(), gm.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(s.second.size(), gm.selected_positions.size());

            std::vector<double> pos;
            for (auto &&i : s.first)
                pos.push_back(i.first);
            // Make sure all positions come out ok
            for (std::size_t i = 0; i < pos.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(pos[i], m.neutral_positions[i]);
                    BOOST_REQUIRE_EQUAL(pos[i], gm.neutral_positions[i]);
                }
            pos.clear();
            for (auto &&i : s.second)
                pos.push_back(i.first);
            for (std::size_t i = 0; i < pos.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(pos[i], m.selected_positions[i]);
                    BOOST_REQUIRE_EQUAL(pos[i], gm.selected_positions[i]);
                }
            auto sums = fwdpp::row_sums(m);
            auto gm_sums = fwdpp::row_sums(gm);
            BOOST_REQUIRE_EQUAL(sums.first.size(), m.neutral_positions.size());
            BOOST_REQUIRE_EQUAL(sums.second.size(),
                                m.selected_positions.size());
            // Check that haplotype and genotype matrices
            // have same row sums in same order
            for (std::size_t i = 0; i < sums.first.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(sums.first[i], gm_sums.first[i]);
                }
            for (std::size_t i = 0; i < sums.second.size(); ++i)
                {
                    BOOST_REQUIRE_EQUAL(sums.second[i], gm_sums.second[i]);
                }
            // Make sure row sums check out for neutral sites.
            // Only compare to haplotype matrix. If we've made it this far,
            // we know m == gm in terms of row sums.
            for (std::size_t c = 0; c < sums.first.size(); ++c)
                {
                    unsigned sum2 = std::count(s.first[c].second.begin(),
                                               s.first[c].second.end(), '1');
                    BOOST_REQUIRE_EQUAL(sums.first[c], sum2);
                }
            // Make sure column sums check out for selected sites
            for (std::size_t c = 0; c < sums.second.size(); ++c)
                {
                    unsigned sum2 = std::count(s.second[c].second.begin(),
                                               s.second[c].second.end(), '1');
                    BOOST_REQUIRE_EQUAL(sums.second[c], sum2);
                }
            // Now, row sums.
            sums = fwdpp::row_sums(m);
            gm_sums = fwdpp::row_sums(gm);
            // Compare row sums from haplotype and genotype matrix.
            for (std::size_t hi = 0, gi = 0; hi < sums.first.size();
                 hi += 1, gi += 1)
                {
                    // Neutral sites
                    BOOST_REQUIRE_EQUAL(sums.first[hi], gm_sums.first[gi]);
                }
            for (std::size_t hi = 0, gi = 0; hi < sums.second.size();
                 hi += 1, gi += 1)
                {
                    // Selected sites
                    BOOST_REQUIRE_EQUAL(sums.second[hi], gm_sums.second[gi]);
                }
            // Check neutral sites in haplotype matrix to the sample
            for (std::size_t r = 0; r < sums.first.size(); ++r)
                {
                    unsigned sum2 = std::count(s.first[r].second.begin(),
                                               s.first[r].second.end(), '1');
                    BOOST_REQUIRE_EQUAL(sums.first[r], sum2);
                }
            // Check selected sites
            for (std::size_t r = 0; r < sums.second.size(); ++r)
                {
                    unsigned sum2 = std::count(s.second[r].second.begin(),
                                               s.second[r].second.end(), '1');
                    BOOST_REQUIRE_EQUAL(sums.second[r], sum2);
                }
            // Evolve pop 100 more generations
            generation = simulate_slocuspop(pop, rng, generation, 100);
        }
}

BOOST_AUTO_TEST_CASE(multilocus_matrix_test)
{
    mlocuspop_popgenmut_fixture mpf;
    simulate_mlocuspop(mpf.pop, mpf.rng, mpf.mutmodels, mpf.recmodels,
                       mlocuspop_popgenmut_fixture::multilocus_additive(),
                       mpf.mu, mpf.rbw, mpf.generation, 10000);
    std::vector<std::size_t> indlist;
    // Sample a LOT of individuals
    for (std::size_t i = 100; i < 750; i += 5)
        indlist.push_back(i);
    for (auto t = 0; t < 100; ++t)
        {
            auto keys = mutation_keys(mpf.pop, indlist, true, true);
            std::sort(keys.first.begin(), keys.first.end(),
                      [&mpf](const std::pair<std::size_t, fwdpp::uint_t> &a,
                             const std::pair<std::size_t, fwdpp::uint_t> &b) {
                          return mpf.pop.mutations[a.first].pos
                                 < mpf.pop.mutations[b.first].pos;
                      });
            std::sort(keys.second.begin(), keys.second.end(),
                      [&mpf](const std::pair<std::size_t, fwdpp::uint_t> &a,
                             const std::pair<std::size_t, fwdpp::uint_t> &b) {
                          return mpf.pop.mutations[a.first].pos
                                 < mpf.pop.mutations[b.first].pos;
                      });
            // In order to be comparable to a sample taken from the population,
            // we must remove sites fixed in these inidivudals
            keys.first.erase(
                std::remove_if(
                    keys.first.begin(), keys.first.end(),
                    [&indlist](
                        const std::pair<std::size_t, fwdpp::uint_t> &a) {
                        return a.second == 2 * indlist.size();
                    }),
                keys.first.end());
            keys.second.erase(
                std::remove_if(
                    keys.second.begin(), keys.second.end(),
                    [&indlist](
                        const std::pair<std::size_t, fwdpp::uint_t> &a) {
                        return a.second == 2 * indlist.size();
                    }),
                keys.second.end());
            auto m
                = haplotype_matrix(mpf.pop, indlist, keys.first, keys.second);
            auto gm
                = genotype_matrix(mpf.pop, indlist, keys.first, keys.second);

            auto sums = fwdpp::row_sums(m);
            auto gm_sums = fwdpp::row_sums(gm);
            // Compare row sums from haplotype and genotype matrix.
            for (std::size_t hi = 0, gi = 0; hi < sums.first.size();
                 hi += 1, gi += 1)
                {
                    // Neutral sites
                    BOOST_REQUIRE_EQUAL(sums.first[hi], gm_sums.first[gi]);
                }
            for (std::size_t hi = 0, gi = 0; hi < sums.second.size();
                 hi += 1, gi += 1)
                {
                    // Selected sites
                    BOOST_REQUIRE_EQUAL(sums.second[hi], gm_sums.second[gi]);
                }

            // Get a sample
            auto s = fwdpp::sample_separate(mpf.pop, indlist, true);
            // multi-locus samples are tricky, so let's simplify
            fwdpp::sample_t neutral, selected;
            std::size_t nlen = 0, slen = 0;
            for (auto &&si : s)
                {
                    nlen += si.first.size();
                    slen += si.second.size();
                    neutral.insert(neutral.end(), si.first.begin(),
                                   si.first.end());
                    selected.insert(selected.end(), si.second.begin(),
                                    si.second.end());
                }
            BOOST_REQUIRE_EQUAL(nlen, m.neutral.size() / m.ncol);
            BOOST_REQUIRE_EQUAL(slen, m.selected.size() / m.ncol);
            BOOST_REQUIRE_EQUAL(sums.first.size(), neutral.size());
            BOOST_REQUIRE_EQUAL(sums.second.size(), selected.size());
            // Check neutral sites in haplotype matrix to the sample
            for (std::size_t r = 0; r < sums.first.size(); ++r)
                {
                    unsigned sum2 = std::count(neutral[r].second.begin(),
                                               neutral[r].second.end(), '1');
                    BOOST_REQUIRE_EQUAL(sums.first[r], sum2);
                }
            // Check selected sites
            for (std::size_t r = 0; r < sums.second.size(); ++r)
                {
                    unsigned sum2 = std::count(selected[r].second.begin(),
                                               selected[r].second.end(), '1');
                    BOOST_REQUIRE_EQUAL(sums.second[r], sum2);
                }
            simulate_mlocuspop(
                mpf.pop, mpf.rng, mpf.mutmodels, mpf.recmodels,
                mlocuspop_popgenmut_fixture::multilocus_additive(), mpf.mu,
                mpf.rbw, mpf.generation, 100);
        }
}
