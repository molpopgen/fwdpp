/*!
  \file siteDepFitness.cc
  \brief Unit tests of KTfwd::site_dependent_fitness
  \ingroup unit
*/

#include <config.h>
#include <cmath>
#include <boost/test/unit_test.hpp>
#include <testsuite/fixtures/fwdpp_fixtures.hpp>
#include <testsuite/util/custom_dip.hpp>
#include <fwdpp/fitness_models.hpp>
#include <memory>
using mut = KTfwd::mutation;

BOOST_FIXTURE_TEST_SUITE(test_site_dependent_fitness,
                         standard_empty_single_deme_fixture)
/*
  This test creates a situation
  where gamete1 has a selected mutation and gamete2 does not.
  If issue #8 were to cause fitness to be mis-calculated,
  then this test will fail.

  However, it does not.  Even with the bug, the remaining bit of the function
  gets the calculation right.  Yay!
*/
BOOST_AUTO_TEST_CASE(simple_multiplicative1)
{
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    g1.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);

    gcont_t g{ g1, g2 };
    double w = KTfwd::site_dependent_fitness()(
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
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, -0.1, 0.5);
    g1.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);

    gcont_t g{ g1, g2 };
    double w = KTfwd::multiplicative_diploid(KTfwd::mtrait())(g[0], g[1],
                                                              mutations);

    BOOST_CHECK_CLOSE(w, -0.05, 1e-8);
}

BOOST_AUTO_TEST_CASE(gss_multiplicative_trait) // Gaussian stab sel
{
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, -0.1, 0.5);
    g1.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);

    gcont_t g{ g1, g2 };
    double optimum = 0.1;
    // Variance in gaussian fitness fxn:
    double VS = 1.0;

    auto gss_closure = [optimum, VS](const double d) {
        // We use d-1 as the genetic value b/c the model
        // is multiplicative, and so we start from a value
        // of 1.0.  Thus, to center the trait values on 0.0,
        // subtract 1.0.
        return std::exp(-std::pow((d - 1.0) - optimum, 2.0) / (2.0 * VS));
    };

    // the above closure changes the behavior
    // so that the genetic value is a trait value
    // mapped to fitness via a "gss" model with
    // params optimum and VS
    double w
        = KTfwd::multiplicative_diploid(gss_closure)(g[0], g[1], mutations);

    BOOST_CHECK_CLOSE(w, std::exp(-std::pow(-0.05 - optimum, 2.0) / (2. * VS)),
                      1e-8);
}
/*
  g2 has it, g1 does not
*/
BOOST_AUTO_TEST_CASE(simple_multiplicative2)
{
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    g2.smutations.emplace_back(0);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 0);
    BOOST_CHECK_EQUAL(g2.smutations.size(), 1);

    gcont_t g{ g1, g2 };

    double w = KTfwd::site_dependent_fitness()(
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
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    g1.smutations.emplace_back(0);
    g2.smutations.emplace_back(0);

    BOOST_CHECK_EQUAL(g1.smutations.size(), 1);
    BOOST_CHECK_EQUAL(g2.smutations.size(), 1);

    gcont_t g{ g1, g2 };

    double w = KTfwd::site_dependent_fitness()(
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
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=0.5 (but we won't use
    // the dominance...)
    mutations.emplace_back(0.1, 0.1, 0.5);
    mutations.emplace_back(0.2, 0.1, 0.5);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    double w = KTfwd::site_dependent_fitness()(
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
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=1.0
    mutations.emplace_back(0.1, 0.1, 1);
    mutations.emplace_back(0.2, 0.1, 1);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    double w = KTfwd::additive_diploid()(g[0], g[1], mutations);
    BOOST_CHECK_EQUAL(w, 1.2);
}

BOOST_AUTO_TEST_CASE(simple_additive_trait)
{
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=1.0
    mutations.emplace_back(0.1, 0.1, 1);
    // s=-0.2 here
    mutations.emplace_back(0.2, -0.2, 1);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    double w = KTfwd::additive_diploid(KTfwd::atrait())(g[0], g[1], mutations);
    BOOST_CHECK_EQUAL(w, -0.1);
}

BOOST_AUTO_TEST_CASE(stateful_additive_trait)
{
    KTfwd::gamete g1(1), g2(1);

    // add mutation at position 0.1, s=0.1,n=1,dominance=1.0
    mutations.emplace_back(0.1, 0.1, 1);
    // s=-0.2 here
    mutations.emplace_back(0.2, -0.2, 1);
    g1.smutations.emplace_back(0);
    g1.smutations.emplace_back(1);
    BOOST_CHECK_EQUAL(g1.smutations.size(), 2);

    gcont_t g{ g1, g2 };

    // There is no logic to this model.
    // More reasonable use cases
    // Would be time-dependent change in
    // optima, etc.
    std::shared_ptr<int> i(new int(0));
    struct stateful_mapping_to_fitness
    {
        // Stateful stuff is harder
        // b/c the object is move-constructed
        // into the genetic value calculator.
        // Thus, we need to hold pointers
        // to external data.
		// Here, we use std::shared_ptr
		// because we do things the right way
		// and don't use bare pointers.
		// A shared_ptr/weak_ptr pairing
		// would also be cool.
        const std::shared_ptr<const int> i;
        stateful_mapping_to_fitness(std::shared_ptr<const int> i_) : i(i_) {}
        inline double
        operator()(const double) const
        {
            if (*(this->i))
                return 1.0;
            return 0.0;
        }
    };

    stateful_mapping_to_fitness sm(i);
	//Here, a COPY of sm is move-constructed
	//into an instance of additive_diploid,
	//which is why we need to use a pointer
	//for our stateful object above.
    auto a = KTfwd::additive_diploid(sm);
    double w = a(g[0], g[1], mutations);
    BOOST_CHECK_EQUAL(w, 0.0);

    *i = 1;
    w = a(g[0], g[1], mutations);
    BOOST_CHECK_EQUAL(w, 1.0);
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
        = KTfwd::traits::fitness_fxn_t<dipvector_t, gcont_t, mcont_t>;

    static_assert(!std::is_same<void, fitness_model_t>::value,
                  "Fitness function signature evaluated to void.");

    // Do bare minimum setup to be able to make calls to functions
    gametes.emplace_back(200);
    mutations.emplace_back(1.0, 0.1); // position 1.0, s = 0.1
    gametes[0].smutations.push_back(0);

    {
        // Test reassignment of the SAME fitness model type
        // Multiplicative model first

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit
            = std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3);

        auto a = dipfit(dipvector_t::value_type{ 0, 0 }, gametes, mutations);
        // Now, reassign it with scaling = 2.
        dipfit
            = std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3, 2.);
        auto b = dipfit(dipvector_t::value_type{ 0, 0 }, gametes, mutations);

        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Test reassignment of the SAME fitness model type
        // Now the additive model

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit
            = std::bind(KTfwd::additive_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3);
        auto a = dipfit(dipvector_t::value_type{ 0, 0 }, gametes, mutations);
        // Now, reassign it with scaling = 2.
        dipfit = std::bind(KTfwd::additive_diploid(), std::placeholders::_1,
                           std::placeholders::_2, std::placeholders::_3, 2.);
        auto b = dipfit(dipvector_t::value_type{ 0, 0 }, gametes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Convert a multiplicative model to an additive model
        // The application of this is a program using the library
        // that wants to decide which model to use based on parameters
        // passed in by a user.

        fitness_model_t dipfit
            = std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3);
        auto a = dipfit(dipvector_t::value_type{ 0, 0 }, gametes, mutations);
        // Now, reassign it to addtive model with scaling = 2.
        dipfit = std::bind(KTfwd::additive_diploid(), std::placeholders::_1,
                           std::placeholders::_2, std::placeholders::_3, 2.);
        auto b = dipfit(dipvector_t::value_type{ 0, 0 }, gametes, mutations);
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
        = KTfwd::traits::fitness_fxn_t<std::vector<custom_diploid_testing_t>,
                                       gcont_t, mcont_t>;

    static_assert(!std::is_same<void, fitness_model_t>::value,
                  "Fitness function signature evaluated to void.");

    // Bare minimum setup for testing
    std::vector<custom_diploid_testing_t> cdiploids(
        1, custom_diploid_testing_t(0, 0));
    gametes.emplace_back(200);
    mutations.emplace_back(1.0, 0.1); // position 1.0, s = 0.1
    gametes[0].smutations.push_back(0);

    {
        // Test reassignment of the SAME fitness model type
        // Multiplicative model first

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit
            = std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3);
        auto a = dipfit(cdiploids[0], gametes, mutations);

        // Now, reassign it with scaling = 2.
        dipfit
            = std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3, 2.);
        auto b = dipfit(cdiploids[0], gametes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Test reassignment of the SAME fitness model type
        // Now the additive model

        // assign a fitness model with default scaling = 1.
        fitness_model_t dipfit
            = std::bind(KTfwd::additive_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3);
        auto a = dipfit(cdiploids[0], gametes, mutations);
        // Now, reassign it with scaling = 2.
        dipfit = std::bind(KTfwd::additive_diploid(), std::placeholders::_1,
                           std::placeholders::_2, std::placeholders::_3, 2.);
        auto b = dipfit(cdiploids[0], gametes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }

    {
        // Convert a multiplicative model to an additive model
        // The application of this is a program using the library
        // that wants to decide which model to use based on parameters
        // passed in by a user.

        fitness_model_t dipfit
            = std::bind(KTfwd::multiplicative_diploid(), std::placeholders::_1,
                        std::placeholders::_2, std::placeholders::_3);
        auto a = dipfit(cdiploids[0], gametes, mutations);
        // Now, reassign it to addtive model with scaling = 2.
        dipfit = std::bind(KTfwd::additive_diploid(), std::placeholders::_1,
                           std::placeholders::_2, std::placeholders::_3, 2.);
        auto b = dipfit(cdiploids[0], gametes, mutations);
        BOOST_CHECK_CLOSE(a - 1.0, 0.5 * (b - 1.0), 1e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END()
