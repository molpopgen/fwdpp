/*!
  \file mlocusCrossoverTest.cc
  \ingroup unit
  \brief Tests KTfwd::fwdpp_internal::multilocus_rec
*/
#include <config.h>
#include <iostream>
// For this unit test, this symbol eliminates the mutation-related part of
// KTfwd::fwdpp_internal::multiloc_rec_mut,
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
    auto x = KTfwd::make_poisson_interlocus_rec(r, rates.data(), rates.size());
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
    // KTfwd::fwdpp_internal::gamete_lookup_table(gametes,mutations);
    auto mutation_recycling_bin
        = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
    auto gamete_recycling_bin
        = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
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

    auto fake_mut_pol
        = [](std::queue<std::size_t> &, decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                          decltype(mutations) &)>>
        mutation_models(3, fake_mut_pol);

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = KTfwd::fwdpp_internal::multilocus_rec_mut(
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
                                            std::vector<double>{ MVAL },
                                            std::vector<double>{ MVAL } };

    // We use these to "fake" what we want to happen between loci.
    std::vector<double> r_bw_loci = { 1., 0. };
    std::vector<diploid_t> diploids({ diploid });

    // auto gamete_lookup =
    // KTfwd::fwdpp_internal::gamete_lookup_table(gametes,mutations);
    auto mutation_recycling_bin
        = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
    auto gamete_recycling_bin
        = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
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

    auto fake_mut_pol
        = [](std::queue<std::size_t> &, decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                          decltype(mutations) &)>>
        mutation_models{ fake_mut_pol, fake_mut_pol };

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = KTfwd::fwdpp_internal::multilocus_rec_mut(
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
        std::vector<double>{ 0.1, 0.45, MVAL }, std::vector<double>{ MVAL },
        std::vector<double>{ MVAL }
    };

    // We use these to "fake" what we want to happen between loci.
    std::vector<double> r_bw_loci = { 1., 1. };
    std::vector<diploid_t> diploids({ diploid });

    // auto gamete_lookup =
    // KTfwd::fwdpp_internal::gamete_lookup_table(gametes,mutations);
    auto mutation_recycling_bin
        = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
    auto gamete_recycling_bin
        = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
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
    auto fake_mut_pol
        = [](std::queue<std::size_t> &, decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                          decltype(mutations) &)>>
        mutation_models{ fake_mut_pol, fake_mut_pol };

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = KTfwd::fwdpp_internal::multilocus_rec_mut(
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

    auto mutation_recycling_bin
        = KTfwd::fwdpp_internal::make_mut_queue(mcounts);
    auto gamete_recycling_bin
        = KTfwd::fwdpp_internal::make_gamete_queue(gametes);
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
    auto fake_mut_pol
        = [](std::queue<std::size_t> &, decltype(mutations) &) { return 0; };
    std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                          decltype(mutations) &)>>
        mutation_models(3, fake_mut_pol);

    double mu[3] = { 0.0, 0.0, 0.0 };

    auto offspring = KTfwd::fwdpp_internal::multilocus_rec_mut(
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
