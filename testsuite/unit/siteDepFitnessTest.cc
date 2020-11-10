/*!
  \file siteDepFitnessTest.cc
  \brief Unit tests of fwdpp::site_dependent_fitness
  \ingroup unit
*/

#include <config.h>
#include <cmath>
#include <boost/test/unit_test.hpp>
#include <testsuite/fixtures/fwdpp_fixtures.hpp>
#include <testsuite/util/custom_dip.hpp>
#include <fwdpp/fitness_models.hpp>
#include <memory>
#include "../../examples/custom_mutation_example.hpp"

using mut = mutation;

BOOST_FIXTURE_TEST_SUITE(test_site_dependent_fitness,
                         standard_empty_single_deme_fixture)
/*
  This test creates a situation
  where haploid_genome1 has a selected mutation and haploid_genome2 does not.
  If issue #8 were to cause fitness to be mis-calculated,
  then this test will fail.

  However, it does not.  Even with the bug, the remaining bit of the function
  gets the calculation right.  Yay!
*/
BOOST_AUTO_TEST_CASE(simple_multiplicative1)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    g1.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);

    gcont_t g{ g1, g2 };
    double w = fwdpp::site_dependent_fitness()(
        g[0], g[1], mutations,
        [&](double &fitness, const mut &__mut) {
            fitness *= std::pow(1. + __mut.s, 2.);
        },
        [&](double &fitness, const mut &__mut) { fitness *= (1. + __mut.s); },
        1.);

    BOOST_CHECK_EQUAL(w, 1.1);
}

BOOST_AUTO_TEST_CASE(simple_multiplicative_trait)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, -0.1, 0.5);
    g1.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);

    gcont_t g{ g1, g2 };
    double w = fwdpp::multiplicative_diploid(fwdpp::trait(1.))(g[0], g[1],
                                                               mutations);

    BOOST_CHECK_CLOSE(w, -0.05, 1e-8);
}

BOOST_AUTO_TEST_CASE(simple_multiplicative_trait_custom_fxn)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, -0.1, 0.5);
    g1.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);

    gcont_t g{ g1, g2 };
    double w = fwdpp::multiplicative_diploid(fwdpp::trait(1.),[](double d) { return d-1.; })(g[0], g[1],
                                                               mutations);

    BOOST_CHECK_CLOSE(w, -0.05, 1e-8);
}
/*
  g2 has it, g1 does not
*/
BOOST_AUTO_TEST_CASE(simple_multiplicative2)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    g2.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 0);
    BOOST_CHECK_EQUAL(g2.smutations.size(), 1);

    gcont_t g{ g1, g2 };

    double w = fwdpp::site_dependent_fitness()(
        g[0], g[1], mutations,
        [&](double &fitness, const mut &__mut) {
            fitness *= std::pow(1. + __mut.s, 2.);
        },
        [&](double &fitness, const mut &__mut) { fitness *= (1. + __mut.s); },
        1.);
    BOOST_CHECK_EQUAL(w, 1.1);
}

/*
  Both have it
*/
BOOST_AUTO_TEST_CASE(simple_multiplicative3)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    g1.smutations.emplace_back(0);
    g2.smutations.emplace_back(0);

    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);
    BOOST_CHECK_EQUAL(g2.smutations.size(), 1);

    gcont_t g{ g1, g2 };

    double w = fwdpp::site_dependent_fitness()(
        g[0], g[1], mutations,
        [&](double &fitness, const mut &__mut) {
            fitness *= std::pow(1. + __mut.s, 2.);
        },
        [&](double &fitness, const mut &__mut) { fitness *= (1. + __mut.s); },
        1.);
    BOOST_CHECK_EQUAL(w, 1.1 * 1.1);
}

/*
  Now, g1 has 2 mutations
*/
BOOST_AUTO_TEST_CASE(simple_multiplicative4)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    mutations.emplace_back(0.2, 0.1, 0.5);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    double w = fwdpp::site_dependent_fitness()(
        g[0], g[1], mutations,
        [&](double &fitness, const mut &__mut) {
            fitness *= std::pow(1. + __mut.s, 2.);
        },
        [&](double &fitness, const mut &__mut) { fitness *= (1. + __mut.s); },
        1.);
    BOOST_CHECK_EQUAL(w, 1.1 * 1.1);
}

BOOST_AUTO_TEST_CASE(simple_additive_1)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=1.0
    mutations.emplace_back(0.1, 0.1, 1);
    mutations.emplace_back(0.2, 0.1, 1);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    auto wfxn = fwdpp::additive_diploid(fwdpp::fitness(1.));
    BOOST_REQUIRE_EQUAL(wfxn.gvalue_is_trait, false);
    BOOST_REQUIRE_EQUAL(wfxn.gvalue_is_fitness, true);
    double w = wfxn(g[0], g[1], mutations);
    BOOST_CHECK_EQUAL(w, 1.2);
}

BOOST_AUTO_TEST_CASE(simple_additive_trait)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=1.0
    mutations.emplace_back(0.1, 0.1, 1);
    // s=-0.2 here
    mutations.emplace_back(0.2, -0.2, 1);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    auto wfxn = fwdpp::additive_diploid(fwdpp::trait(1.));
    BOOST_REQUIRE_EQUAL(wfxn.gvalue_is_trait, true);
    BOOST_REQUIRE_EQUAL(wfxn.gvalue_is_fitness, false);
    auto w = wfxn(g[0], g[1], mutations);
    BOOST_CHECK_EQUAL(w, -0.1);
}

BOOST_AUTO_TEST_CASE(simple_additive_trait_custom_fxn)
{
    fwdpp::haploid_genome g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=1.0
    mutations.emplace_back(0.1, 0.1, 1);
    // s=-0.2 here
    mutations.emplace_back(0.2, -0.2, 1);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    auto wfxn = fwdpp::additive_diploid(fwdpp::trait(1.),[](double d){return d;});
    BOOST_REQUIRE_EQUAL(wfxn.gvalue_is_trait, true);
    BOOST_REQUIRE_EQUAL(wfxn.gvalue_is_fitness, false);
    auto w = wfxn(g[0], g[1], mutations);
    BOOST_CHECK_EQUAL(w, -0.1);
}

/*
  API checks on fitness policies.

  Below, we test the ability to assign bound fitness models to variables
  and then reassign them.

  We make use of types defined by population objects from the "sugar" layer.
*/

BOOST_AUTO_TEST_CASE(reassign_test_1)
{
    // This is the expected type of a fitness policy for non-custom diploids.
    // We use fwdpp's type_traits header to pull this out.
    using fitness_model_t
        = fwdpp::traits::fitness_fxn_t<dipvector_t, gcont_t, mcont_t>;

    static_assert(!std::is_same<void, fitness_model_t>::value,
                  "Fitness function signature evaluated to void.");

    // Do bare minimum setup to be able to make calls to functions
    haploid_genomes.emplace_back(200);
    mutations.emplace_back(1.0, 0.1); // position 1.0, s = 0.1
    haploid_genomes[0].smutations.push_back(0);

    {
        // Test reassignment of the SAME fitness model type
        // Multiplicative model first

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit
            = std::bind(fwdpp::multiplicative_diploid(fwdpp::fitness(1.)),
                        std::placeholders::_1, std::placeholders::_2,
                        std::placeholders::_3);

        auto a = dipfit(dipvector_t::value_type{ 0, 0 }, haploid_genomes, mutations);
        // Now, reassign it with scaling = 2.
        dipfit = fwdpp::multiplicative_diploid(fwdpp::fitness(2.));
        auto b = dipfit(dipvector_t::value_type{ 0, 0 }, haploid_genomes, mutations);

        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Test reassignment of the SAME fitness model type
        // Now the additive model

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit = fwdpp::additive_diploid(fwdpp::fitness(1.));
        auto a = dipfit(dipvector_t::value_type{ 0, 0 }, haploid_genomes, mutations);
        // Now, reassign it with scaling = 2.
        dipfit = fwdpp::additive_diploid(fwdpp::fitness(2.));
        auto b = dipfit(dipvector_t::value_type{ 0, 0 }, haploid_genomes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Convert a multiplicative model to an additive model
        // The application of this is a program using the library
        // that wants to decide which model to use based on parameters
        // passed in by a user.

        fitness_model_t dipfit
            = fwdpp::multiplicative_diploid(fwdpp::fitness(1.));
        auto a = dipfit(dipvector_t::value_type{ 0, 0 }, haploid_genomes, mutations);
        // Now, reassign it to addtive model with scaling = 2.
        dipfit = fwdpp::additive_diploid(fwdpp::fitness(2.0));
        auto b = dipfit(dipvector_t::value_type{ 0, 0 }, haploid_genomes, mutations);
        // With only 1 mutation in play, additive & multiplicative will give
        // same result:
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(reassign_test_2)
{
    // Expected signature for fitness policies involving custom diploids
    // We use fwdpp's type_traits header to pull this out.
    using fitness_model_t
        = fwdpp::traits::fitness_fxn_t<std::vector<custom_diploid_testing_t>,
                                       gcont_t, mcont_t>;

    static_assert(!std::is_same<void, fitness_model_t>::value,
                  "Fitness function signature evaluated to void.");

    // Bare minimum setup for testing
    std::vector<custom_diploid_testing_t> cdiploids(
        1, custom_diploid_testing_t(0, 0));
    haploid_genomes.emplace_back(200);
    mutations.emplace_back(1.0, 0.1); // position 1.0, s = 0.1
    haploid_genomes[0].smutations.push_back(0);

    {
        // Test reassignment of the SAME fitness model type
        // Multiplicative model first

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit
            = fwdpp::multiplicative_diploid(fwdpp::fitness(1.));
        auto a = dipfit(cdiploids[0], haploid_genomes, mutations);

        // Now, reassign it with scaling = 2.
        dipfit = fwdpp::multiplicative_diploid(fwdpp::fitness(2.));
        auto b = dipfit(cdiploids[0], haploid_genomes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Test reassignment of the SAME fitness model type
        // Now the additive model

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit = fwdpp::additive_diploid(fwdpp::fitness(1.));
        auto a = dipfit(cdiploids[0], haploid_genomes, mutations);
        // Now, reassign it with scaling = 2.
        dipfit = fwdpp::additive_diploid(fwdpp::fitness(2.));
        auto b = dipfit(cdiploids[0], haploid_genomes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Convert a multiplicative model to an additive model
        // The application of this is a program using the library
        // that wants to decide which model to use based on parameters
        // passed in by a user.

        fitness_model_t dipfit
            = fwdpp::multiplicative_diploid(fwdpp::fitness(1.));
        auto a = dipfit(cdiploids[0], haploid_genomes, mutations);
        // Now, reassign it to addtive model with scaling = 2.
        dipfit = fwdpp::additive_diploid(fwdpp::fitness(2.));
        auto b = dipfit(cdiploids[0], haploid_genomes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END()
