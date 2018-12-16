/*!
  \file mlocusCrossoverTest.cc
  \ingroup unit
  \brief Tests fwdpp::fwdpp_internal::multilocus_rec
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
#include "../fixtures/fwdpp_fixtures.hpp"
#include "../fixtures/multilocus_fixture_deterministic.hpp"

/*
  Initiate random number generation system --
  these tests will not be random, but these objects
  are req'd for function calls to fwdpp's internal
  stuff
*/
gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlxs2);

using diploid_t = dipvector_t;

// /*
//   Setup fxn for 3-locus scenario

//   Set up the following config:
//   g1,l1 = 0.5
//   g1,l2 = 0.75,
//   g2,l1 = 0.25,
//   g2,l2 = 0.9
//   g1,l3 = 1.25, 1.5
//   g2,l3 = 1.1
// */
void
setup3locus2(gcont_t &gametes, mcont_t &mutations, diploid_t &diploid)
{
    // gametes = gcont_t(6,gcont_t::value_type(1));
    gametes.resize(6, gcont_t::value_type(1));
    mutations.clear();

    // To set this up, let's add the mutations:
    mutations.emplace_back(0.5, 0.);
    mutations.emplace_back(0.75, 0.);
    mutations.emplace_back(0.25, 0.);
    mutations.emplace_back(0.9, 0.);
    mutations.emplace_back(1.25, 0.1);
    mutations.emplace_back(1.5, 0.1);
    mutations.emplace_back(1.1, 0.1);

    // put the mutations into gametes
    gametes[0].mutations.push_back(0); // gamete 1, locus 1
    gametes[1].mutations.push_back(2); // gamete 2, locus 1
    gametes[2].mutations.push_back(1);
    gametes[3].mutations.push_back(3);
    gametes[4].mutations.push_back(4);
    gametes[4].mutations.push_back(5);
    gametes[5].mutations.push_back(6);

    // Now, make a diploid
    diploid.clear();
    diploid.emplace_back(std::make_pair(0, 1));
    diploid.emplace_back(std::make_pair(2, 3));
    diploid.emplace_back(std::make_pair(4, 5));
    return;
}

struct always_recombine
{
    inline unsigned
    operator()() const
    {
        return 1u;
    }
};

BOOST_FIXTURE_TEST_SUITE(multilocus_rec_mutTest,
                         standard_empty_multiloc_fixture)

BOOST_AUTO_TEST_CASE(test_flexibility)
{
    std::vector<double> rates{ 1e-3, 0.1 };
    auto x = fwdpp::make_poisson_interlocus_rec(r, rates.data(), rates.size());
    x.emplace_back(always_recombine());
}

BOOST_AUTO_TEST_CASE(three_locus_test_1)
{
    diploid_t diploid; // parent 1
    setup3locus2(gametes, mutations, diploid);
    std::vector<unsigned> mcounts(mutations.size(), 1);
    diploid_t diploid2(diploid); // parent 2

    // positions of x-overs within loci
    auto MVAL = std::numeric_limits<double>::max();
    std::vector<std::vector<double>> rec1{ std::vector<double>{ 0.3, MVAL },
                                           std::vector<double>{ 0.55, 0.8,
                                                                MVAL },
                                           std::vector<double>{ 1.3, MVAL } };

    // We use these to "fake" what we want to happen between loci.
    std::vector<double> r_bw_loci = { 1., 0. };
    std::vector<diploid_t> diploids({ diploid });

    // auto gamete_lookup =
    // fwdpp::fwdpp_internal::gamete_lookup_table(gametes,mutations);
    auto mutation_recycling_bin = fwdpp::make_mut_queue(mcounts);
    auto gamete_recycling_bin = fwdpp::make_gamete_queue(gametes);
    gcont_t::value_type::mutation_container neutral,
        selected; // req'd as of 0.3.3

    std::vector<std::function<std::vector<double>(const gcont_t::value_type &,
                                                  const gcont_t::value_type &,
                                                  const mcont_t &)>>
        recpols{
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[0]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[1]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[2]; }
        };

    std::vector<std::function<unsigned(void)>> interlocus_rec{
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[0]); },
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[1]); }
    };

    auto fake_mut_pol = [](fwdpp::flagged_mutation_queue &,
                           decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          decltype(mutations) &)>>
        mutation_models(3, fake_mut_pol);

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        r, diploid, diploid2, mutation_recycling_bin, gamete_recycling_bin,
        recpols, interlocus_rec, 0, 0, gametes, mutations, neutral, selected,
        &mu[0], mutation_models);

    BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(), 0);
    BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(), 0);
    BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(), 1);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[2].first].mutations[0]].pos,
                      1.25);
}

BOOST_AUTO_TEST_CASE(three_locus_test_2)
{
    diploid_t diploid; // parent 1
    setup3locus2(gametes, mutations, diploid);
    std::vector<unsigned> mcounts(mutations.size(), 1);
    diploid_t diploid2(diploid); // parent 2

    // positions of x-overs within loci
    auto MVAL = std::numeric_limits<double>::max();
    auto rec1
        = std::vector<std::vector<double>>{ std::vector<double>{ 0.45, MVAL },
                                            std::vector<double>{},
                                            std::vector<double>{} };

    // We use these to "fake" what we want to happen between loci.
    std::vector<double> r_bw_loci = { 1., 0. };
    std::vector<diploid_t> diploids({ diploid });

    // auto gamete_lookup =
    // fwdpp::fwdpp_internal::gamete_lookup_table(gametes,mutations);
    auto mutation_recycling_bin = fwdpp::make_mut_queue(mcounts);
    auto gamete_recycling_bin = fwdpp::make_gamete_queue(gametes);
    gcont_t::value_type::mutation_container neutral,
        selected; // req'd as of 0.3.3

    std::vector<std::function<std::vector<double>(const gcont_t::value_type &,
                                                  const gcont_t::value_type &,
                                                  const mcont_t &)>>
        recpols{
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[0]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[1]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[2]; }
        };

    std::vector<std::function<unsigned(void)>> interlocus_rec{
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[0]); },
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[1]); }
    };

    auto fake_mut_pol = [](fwdpp::flagged_mutation_queue &,
                           decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          decltype(mutations) &)>>
        mutation_models{ fake_mut_pol, fake_mut_pol };

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        r, diploid, diploid2, mutation_recycling_bin, gamete_recycling_bin,
        recpols, interlocus_rec, 0, 0, gametes, mutations, neutral, selected,
        &mu[0], mutation_models);

    BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(), 0);
    BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(), 1);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[1].first].mutations[0]].pos,
                      0.75);
    BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(), 2);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[2].first].mutations[0]].pos,
                      1.25);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[2].first].mutations[1]].pos,
                      1.5);
}

BOOST_AUTO_TEST_CASE(three_locus_test_3)
{
    diploid_t diploid; // parent 1
    setup3locus2(gametes, mutations, diploid);
    std::vector<unsigned> mcounts(mutations.size(), 1);
    diploid_t diploid2(diploid); // parent 2

    // positions of x-overs within loci
    auto MVAL = std::numeric_limits<double>::max();
    auto rec1 = std::vector<std::vector<double>>{
        std::vector<double>{ 0.1, 0.45, MVAL }, std::vector<double>{},
        std::vector<double>{}
    };

    // We use these to "fake" what we want to happen between loci.
    std::vector<double> r_bw_loci = { 1., 1. };
    std::vector<diploid_t> diploids({ diploid });

    // auto gamete_lookup =
    // fwdpp::fwdpp_internal::gamete_lookup_table(gametes,mutations);
    auto mutation_recycling_bin = fwdpp::make_mut_queue(mcounts);
    auto gamete_recycling_bin = fwdpp::make_gamete_queue(gametes);
    gcont_t::value_type::mutation_container neutral,
        selected; // req'd as of 0.3.3

    std::vector<std::function<std::vector<double>(const gcont_t::value_type &,
                                                  const gcont_t::value_type &,
                                                  const mcont_t &)>>
        recpols{
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[0]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[1]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[2]; }
        };

    std::vector<std::function<unsigned(void)>> interlocus_rec{
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[0]); },
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[1]); }
    };
    auto fake_mut_pol = [](fwdpp::flagged_mutation_queue &,
                           decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          decltype(mutations) &)>>
        mutation_models{ fake_mut_pol, fake_mut_pol };

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        r, diploid, diploid2, mutation_recycling_bin, gamete_recycling_bin,
        recpols, interlocus_rec, 0, 0, gametes, mutations, neutral, selected,
        &mu[0], mutation_models);

    BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(), 2);
    BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(), 1);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[1].first].mutations[0]].pos,
                      0.9);
    BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(), 2);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[2].first].mutations[0]].pos,
                      1.25);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[2].first].mutations[1]].pos,
                      1.5);
}

BOOST_AUTO_TEST_CASE(three_locus_test_4)
{
    diploid_t diploid; // parent 1
    setup3locus2(gametes, mutations, diploid);
    std::vector<unsigned> mcounts(mutations.size(), 1);
    diploid_t diploid2(diploid); // parent 2

    // positions of x-overs within loci
    auto MVAL = std::numeric_limits<double>::max();
    auto rec1 = std::vector<std::vector<double>>{
        std::vector<double>{ 0.1, 0.45, 0.51, MVAL },
        std::vector<double>{ 0.6, 0.7, 0.99, MVAL },
        std::vector<double>{ 1, 1.2, 1.3, 1.55, MVAL }
    };

    // We use these to "fake" what we want to happen between loci.
    std::vector<double> r_bw_loci = { 1., 1. };
    std::vector<diploid_t> diploids({ diploid });

    auto mutation_recycling_bin = fwdpp::make_mut_queue(mcounts);
    auto gamete_recycling_bin = fwdpp::make_gamete_queue(gametes);
    gcont_t::value_type::mutation_container neutral,
        selected; // req'd as of 0.3.3

    std::vector<std::function<std::vector<double>(const gcont_t::value_type &,
                                                  const gcont_t::value_type &,
                                                  const mcont_t &)>>
        recpols{
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[0]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[1]; },
            [&rec1](const gcont_t::value_type &, const gcont_t::value_type &,
                    const mcont_t &) { return rec1[2]; }
        };

    std::vector<std::function<unsigned(void)>> interlocus_rec{
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[0]); },
        [&r_bw_loci]() { return static_cast<unsigned>(r_bw_loci[1]); }
    };
    auto fake_mut_pol = [](fwdpp::flagged_mutation_queue &,
                           decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          decltype(mutations) &)>>
        mutation_models(3, fake_mut_pol);

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        r, diploid, diploid2, mutation_recycling_bin, gamete_recycling_bin,
        recpols, interlocus_rec, 0, 0, gametes, mutations, neutral, selected,
        &mu[0], mutation_models);

    BOOST_CHECK_EQUAL(gametes[offspring[0].first].mutations.size(), 2);
    BOOST_CHECK_EQUAL(gametes[offspring[1].first].mutations.size(), 1);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[1].first].mutations[0]].pos,
                      0.75);
    BOOST_CHECK_EQUAL(gametes[offspring[2].first].mutations.size(), 2);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[2].first].mutations[0]].pos,
                      1.1);
    BOOST_CHECK_EQUAL(mutations[gametes[offspring[2].first].mutations[1]].pos,
                      1.25);
}
BOOST_AUTO_TEST_SUITE_END()

// Below is a different set of tests based
// on a cleaner fixture, which is also used
// in testing of tree sequence recording from
// a multi-locus population.
BOOST_FIXTURE_TEST_SUITE(test_multilocus_recombination,
                         multilocus_fixture_deterministic)

BOOST_AUTO_TEST_CASE(test_transmission)
// The expected genotype of all "first" gametes in offspring
// should be identical to parent1/gamete1 at all loci except
// on the interval [1.5,2) and [3.5,4)
{
    mutate_parent();
    // This loop essentially unit-tests add_mutation (again).
    for (std::size_t i = 0; i < nloci; ++i)
        {
            auto &locus = pop.diploids[0][i];
            BOOST_REQUIRE_EQUAL(pop.gametes[locus.first].mutations.size(), 1);
            BOOST_REQUIRE_EQUAL(
                pop.mutations[pop.gametes[locus.first].mutations[0]].pos,
                pop.locus_boundaries[i].first);
        }
    std::vector<double> mu(nloci, 0.);
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          poptype::mcont_t &)>>
        dummy_mmodels;
    int generation = 0;
    for (unsigned i = 0; i < nloci; ++i)
        {
            dummy_mmodels.push_back([this, &generation,
                                     i](fwdpp::flagged_mutation_queue &recbin,
                                        poptype::mcont_t &mutations) {
                return fwdpp::infsites_popgenmut(
                    recbin, mutations, rng.get(), pop.mut_lookup, generation,
                    0.0,
                    [this, i]() { return gsl_ran_flat(rng.get(), i, i + 1); },
                    []() { return 0.0; }, []() { return 0.0; });
            });
        }

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        rng.get(), pop.diploids[0], pop.diploids[1],
        params_no_swap.mutation_recycling_bin,
        params_no_swap.gamete_recycling_bin,
        params_no_swap.generate_breakpoints,
        params_no_swap.interlocus_recombination, 0, 0, pop.gametes,
        pop.mutations, params_no_swap.neutral, params_no_swap.selected,
        mu.data(), dummy_mmodels);

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

BOOST_AUTO_TEST_CASE(test_transmission_with_extra_variants)
// The expected genotype of all "first" gametes in offspring
// should be identical to parent1/gamete1 at all loci except
// on the interval [1.5,2) and [3.5,4)
{
    mutate_parent2();
    for (std::size_t i = 0; i < nloci; ++i)
        {
            auto &locus = pop.diploids[0][i];
            BOOST_REQUIRE_EQUAL(pop.gametes[locus.first].mutations.size(), 2);
            BOOST_REQUIRE_EQUAL(
                pop.mutations[pop.gametes[locus.first].mutations[0]].pos,
                pop.locus_boundaries[i].first);
            BOOST_CHECK_CLOSE(
                pop.mutations[pop.gametes[locus.first].mutations[1]].pos,
                pop.locus_boundaries[i].first + 0.51, 1e-5);
        }
    std::vector<double> mu(nloci, 0.);
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          poptype::mcont_t &)>>
        dummy_mmodels;
    int generation = 0;
    for (unsigned i = 0; i < nloci; ++i)
        {
            dummy_mmodels.push_back([this, &generation,
                                     i](fwdpp::flagged_mutation_queue &recbin,
                                        poptype::mcont_t &mutations) {
                return fwdpp::infsites_popgenmut(
                    recbin, mutations, rng.get(), pop.mut_lookup, generation,
                    0.0,
                    [this, i]() { return gsl_ran_flat(rng.get(), i, i + 1); },
                    []() { return 0.0; }, []() { return 0.0; });
            });
        }

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        rng.get(), pop.diploids[0], pop.diploids[1],
        params_no_swap.mutation_recycling_bin,
        params_no_swap.gamete_recycling_bin,
        params_no_swap.generate_breakpoints,
        params_no_swap.interlocus_recombination, 0, 0, pop.gametes,
        pop.mutations, params_no_swap.neutral, params_no_swap.selected,
        mu.data(), dummy_mmodels);

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

    //Every other locus must have had its internal variant removed
    for (std::size_t i = 0; i < nloci; ++i)
        {
            if (i % 2 == 0.)
                {
                    BOOST_REQUIRE_EQUAL(
                        pop.gametes[offspring[i].first].mutations.size(), 2);
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

BOOST_AUTO_TEST_CASE(test_transmission_swap_1)
{
    mutate_parent();
    for (std::size_t i = 0; i < nloci; ++i)
        {
            auto &locus = pop.diploids[0][i];
            BOOST_REQUIRE_EQUAL(pop.gametes[locus.first].mutations.size(), 1);
            BOOST_REQUIRE_EQUAL(
                pop.mutations[pop.gametes[locus.first].mutations[0]].pos,
                pop.locus_boundaries[i].first);
        }
    std::vector<double> mu(nloci, 0.);
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          poptype::mcont_t &)>>
        dummy_mmodels;
    int generation = 0;
    for (unsigned i = 0; i < nloci; ++i)
        {
            dummy_mmodels.push_back([this, &generation,
                                     i](fwdpp::flagged_mutation_queue &recbin,
                                        poptype::mcont_t &mutations) {
                return fwdpp::infsites_popgenmut(
                    recbin, mutations, rng.get(), pop.mut_lookup, generation,
                    0.0,
                    [this, i]() { return gsl_ran_flat(rng.get(), i, i + 1); },
                    []() { return 0.0; }, []() { return 0.0; });
            });
        }

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        rng.get(), pop.diploids[0], pop.diploids[1],
        params_no_swap.mutation_recycling_bin,
        params_no_swap.gamete_recycling_bin,
        params_no_swap.generate_breakpoints,
        params_no_swap.interlocus_recombination, 1, 0, pop.gametes,
        pop.mutations, params_no_swap.neutral, params_no_swap.selected,
        mu.data(), dummy_mmodels);

    // Check transmission of mutations into offpring's FIRST gamete
    int locus = 0;
    bool expected_result = false;
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

BOOST_AUTO_TEST_CASE(test_transmission_2)
{
    mutate_parent2();
    std::vector<double> mu(nloci, 0.);
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          poptype::mcont_t &)>>
        dummy_mmodels;
    int generation = 0;
    for (unsigned i = 0; i < nloci; ++i)
        {
            dummy_mmodels.push_back([this, &generation,
                                     i](fwdpp::flagged_mutation_queue &recbin,
                                        poptype::mcont_t &mutations) {
                return fwdpp::infsites_popgenmut(
                    recbin, mutations, rng.get(), pop.mut_lookup, generation,
                    0.0,
                    [this, i]() { return gsl_ran_flat(rng.get(), i, i + 1); },
                    []() { return 0.0; }, []() { return 0.0; });
            });
        }

    auto offspring = fwdpp::fwdpp_internal::multilocus_rec_mut(
        rng.get(), pop.diploids[0], pop.diploids[1],
        params_no_swap2.mutation_recycling_bin,
        params_no_swap2.gamete_recycling_bin,
        params_no_swap2.generate_breakpoints,
        params_no_swap2.interlocus_recombination, 1, 0, pop.gametes,
        pop.mutations, params_no_swap2.neutral, params_no_swap2.selected,
        mu.data(), dummy_mmodels);
}
BOOST_AUTO_TEST_SUITE_END()
