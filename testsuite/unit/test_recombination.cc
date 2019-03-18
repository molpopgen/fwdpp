/*!
  \file test_recombination.cc
  \ingroup unit
*/
#include <config.h>
#include <iostream>
// For this unit test, this symbol eliminates the mutation-related part of
// fwdpp::fwdpp_internal::multiloc_rec_mut,
// which means we don't have to write as much boilerplate code to test the more
// complex logic.
// Plus, mutation stuff is unit-tested elsewhere
// #define FWDPP_UNIT_TESTING
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <vector>
#include <gsl/gsl_rng.h>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>

BOOST_AUTO_TEST_SUITE(test_recombination)

struct simple_mutation : public fwdpp::mutation_base
{
    simple_mutation(double p) : fwdpp::mutation_base(p) {}
};

BOOST_AUTO_TEST_CASE(compare_single_vs_multi_locus_rec)
{
    std::vector<fwdpp::gamete> gametes(2, fwdpp::gamete(1));
    std::vector<simple_mutation> mutations;
    std::vector<fwdpp::uint_t> new_mutations, neutral, selected;
    for (int i = 0; i < 10; ++i)
        {
            mutations.emplace_back(static_cast<double>(i) * 0.5);
            if (i % 2 == 0.0)
                {
                    gametes[0].mutations.push_back(i);
                }
            else
                {
                    gametes[1].mutations.push_back(i);
                }
        }
    std::sort(begin(gametes[0].mutations), end(gametes[0].mutations),
              [&mutations](fwdpp::uint_t a, fwdpp::uint_t b) {
                  return mutations[a].pos < mutations[b].pos;
              });
    std::sort(begin(gametes[1].mutations), end(gametes[1].mutations),
              [&mutations](fwdpp::uint_t a, fwdpp::uint_t b) {
                  return mutations[a].pos < mutations[b].pos;
              });
    auto rbin = fwdpp::empty_gamete_queue();
    std::vector<double> breakpoints{
        0.1, 0.25, 2, 2.25, 3.77, 4, 4.5, std::numeric_limits<double>::max()
    };
    auto offspring
        = fwdpp::mutate_recombine(new_mutations, breakpoints, 0, 1, gametes,
                                  mutations, rbin, neutral, selected);
    //Now, set up the multi-locus version
    //Locus boundaries are 2 and 4
    std::vector<std::vector<std::pair<fwdpp::uint_t, fwdpp::uint_t>>> diploids;
    std::vector<fwdpp::gamete> ml_gametes(7, fwdpp::gamete(1));
    for (int i = 0; i < 10; ++i)
        {
            if (i % 2 == 0.0)
                {
                    if (mutations[i].pos < 2)
                        {
                            ml_gametes[0].mutations.push_back(i);
                        }
                    else if (mutations[i].pos < 4)
                        {
                            ml_gametes[1].mutations.push_back(i);
                        }
                    else
                        {
                            ml_gametes[2].mutations.push_back(i);
                        }
                }
            else
                {
                    if (mutations[i].pos < 2)
                        {
                            ml_gametes[3].mutations.push_back(i);
                        }
                    else if (mutations[i].pos < 4)
                        {
                            ml_gametes[4].mutations.push_back(i);
                        }
                    else
                        {
                            ml_gametes[5].mutations.push_back(i);
                        }
                }
        }
    diploids.resize(2);
    diploids[0].resize(3);
    diploids[1].resize(3);
    ml_gametes[6].n = 6;
    diploids[0][0] = { 0, 3 };
    diploids[0][1] = { 1, 4 };
    diploids[0][2] = { 2, 5 };
    diploids[1][0] = { 6, 6 };
    diploids[1][1] = { 6, 6 };
    diploids[1][2] = { 6, 6 };

    std::vector<std::function<unsigned(void)>> interlocus_rec(
        2, []() { return 1; });
    //0.1, 0.25, 2, 2.25, 4, 4.5, std::numeric_limits<double>::max()
    std::vector<std::function<std::vector<double>(void)>> intralocus_rec;
    intralocus_rec.emplace_back([]() {
        return std::vector<double>{ 0.1, 0.25,
                                    std::numeric_limits<double>::max() };
    });
    intralocus_rec.emplace_back([]() {
        return std::vector<double>{ 2.25, 3.77, std::numeric_limits<double>::max() };
    });
    intralocus_rec.emplace_back([]() {
        return std::vector<double>{ 4.5, std::numeric_limits<double>::max() };
    });
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          std::vector<simple_mutation> &)>>
        mmodels;
    for (int i = 0; i < 3; ++i)
        {
            mmodels.emplace_back([](fwdpp::flagged_mutation_queue &,
                                    std::vector<simple_mutation> &mutations) {
                return mutations.size();
            });
        }
    auto mbin = fwdpp::empty_mutation_queue();
    std::vector<double> mu(3, 0.0);
    gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxs2);
    auto ml_offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        r, diploids[0], diploids[1], mbin, rbin, intralocus_rec,
        interlocus_rec, 0, 0, ml_gametes, mutations, neutral, selected,
        mu.data(), mmodels);
    std::vector<fwdpp::uint_t> pooled;
    for (auto g : ml_offspring)
        {
            for (auto k : ml_gametes[g.first].mutations)
                {
                    pooled.push_back(k);
                }
        }
    gsl_rng_free(r);
    BOOST_REQUIRE_EQUAL(gametes[offspring].mutations == pooled, true);
}

BOOST_AUTO_TEST_SUITE_END()
