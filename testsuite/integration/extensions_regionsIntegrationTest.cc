#include <config.h>
#include <cmath>
#include <boost/test/unit_test.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <fwdpp/extensions/regions.hpp>
#include <fwdpp/sample_diploid.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/interlocus_recombination.hpp>
#include <fwdpp/type_traits.hpp>
#include <limits>
#include "../fixtures/sugar_fixtures.hpp"

using namespace fwdpp;

BOOST_FIXTURE_TEST_SUITE(test_regions, slocuspop_popgenmut_fixture)

// Check that extensions::discrete_mut_model::operator() can be bound
// with placeholders, that the resulting type is a valid
// mutation model, and can be passed to fwdpp::sample_diploid
// BOOST_AUTO_TEST_CASE(discrete_mut_model_test_4)
// {
//     // attempt
//     extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 0.5 }, {},
//     {},
//                                       {}, {});
//     fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2> rng(0u);
//     auto mmodel = std::bind(
//         &extensions::discrete_mut_model::operator()
//             < fwdpp::traits::recycling_bin_t<decltype(pop.mutations)>,
//         decltype(pop.mut_lookup), decltype(pop.mutations) >, &dm,
//         std::placeholders::_1, std::placeholders::_2, rng.get(), 0.001, 0.,
//         &generation, std::ref(pop.mut_lookup));
//     static_assert(traits::is_mutation_model<decltype(mmodel),
//     poptype::mcont_t,
//                                             poptype::gcont_t>::value,
//                   "error: type mutation_model is not a dispatchable mutation
//                   " "model type!");
//
//     auto wbar = fwdpp::sample_diploid(
//         rng.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
//         1000, 0.001, mmodel, fwdpp::poisson_xover(rng.get(), 0.001, 0., 2.),
//         fwdpp::multiplicative_diploid(2.), pop.neutral, pop.selected);
//     if (!std::isfinite(wbar))
//         {
//             throw std::runtime_error("wbar not finite");
//         }
// }

// Test the convenience fxn
// BOOST_AUTO_TEST_CASE(discrete_mut_model_test_5)
// {
//     // attempt
//     extensions::discrete_mut_model dm({ 0, 1 }, { 1, 2 }, { 1, 0.5 }, {},
//     {},
//                                       {}, {});
//     fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2> rng(0u);
//     unsigned generation = 0;
//     auto wbar = fwdpp::sample_diploid(
//         rng.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
//         1000, 0.001, extensions::bind_dmm(dm, pop.mutations, pop.mut_lookup,
//                                     rng.get(), 0.001, 0., &generation),
//         fwdpp::poisson_xover(rng.get(), 0.001, 0., 2.),
//         fwdpp::multiplicative_diploid(2.), pop.neutral, pop.selected);
//     if (!std::isfinite(wbar))
//         {
//             throw std::runtime_error("wbar not finite");
//         }
// }

/*
  Now, test discrete_mut_model's constructor that takes labels,
  a feature introduced in 0.4.9.  The purpose of this is to
  assign to mutation_base::xtra, which allows mutations to be integer-labelled.
*/
// BOOST_AUTO_TEST_CASE(discrete_mut_model_test_6)
// // This is an 'integration' test, I guess...
// {
//     // attempt
//     extensions::discrete_mut_model dm(
//         { 0, 1 },   // starts of 'neutral' regions
//         { 1, 2 },   // ends of 'neutral' regions
//         { 1, 0.5 }, // weights on 'neutral' regions
//         {},         // starts of 'selected' regions
//         {},         // stops of 'selected' regions
//         {},         // weights on 'selected' regions
//         {},         // vector of shmodels
//         { 0,
//           1 }, // labels to put on mutations from each of the 'neutral'
//           regions
//         {} // labels to put on mutations from each of the 'selected' regions
//         );
//
//     // now, evolve the population
//     fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2> rng(0u);
//     auto mmodel = std::bind(
//         &extensions::discrete_mut_model::operator()
//             < fwdpp::traits::recycling_bin_t<decltype(pop.mutations)>,
//         decltype(pop.mut_lookup), decltype(pop.mutations) >, &dm,
//         std::placeholders::_1, std::placeholders::_2, rng.get(), 0.001, 0.,
//         &generation, std::ref(pop.mut_lookup));
//     static_assert(traits::is_mutation_model<decltype(mmodel),
//     poptype::mcont_t,
//                                             poptype::gcont_t>::value,
//                   "error: type mutation_model is not a dispatchable mutation
//                   " "model type!");
//     auto wbar = fwdpp::sample_diploid(
//         rng.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
//         1000, 0.01, // mutation rate--high so that there are lots of
//         mutations to
//         // test below...
//         mmodel, fwdpp::poisson_xover(rng.get(), 0.001, 0., 2.),
//         fwdpp::multiplicative_diploid(2.), pop.neutral, pop.selected);
//     if (!std::isfinite(wbar))
//         {
//             throw std::runtime_error("wbar not finite");
//         }
//     // Check that mutations in certain position intervals have the correct
//     // label
//     for (const auto& m : pop.mutations)
//         {
//             if (m.pos < 1)
//                 {
//                     BOOST_REQUIRE_EQUAL(m.xtra, 0);
//                 }
//             else
//                 {
//                     BOOST_REQUIRE_EQUAL(m.xtra, 1);
//                 }
//         }
// }

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(test_multilocus_regions, mlocuspop_popgenmut_fixture)

BOOST_AUTO_TEST_CASE(test_bind_vec_dmm_drm)
/* Test vectors of mutation/recombination regions
 * The recombination region scheme is same as unit
 * test bind_vec_drm_test in file testsuite/unit/extensions_regionsTest.cc
 *
 * The mutation setup is simpler--just neutral mutations.
 *
 * This test uses the multilocus fixture for the testsuite.
 */
{
	BOOST_REQUIRE_EQUAL(vdmm.size(),bound_mmodels.size());
	BOOST_REQUIRE_EQUAL(vdmm.size(),vdrm.size());
    // create a set of bound callbacks.
    // We use the fixture's mu to imply that
    // mutation rate = recombination rate per region.

    // the mutation rates to neutral/selected variants
    std::vector<double> neutral_mutrates(4, 1e-3), selected_mutrates(4, 0.);

    auto interlocus_rec = fwdpp::make_binomial_interlocus_rec(
        rng.get(), rbw.data(), rbw.size());
    double wbar = sample_diploid(
        rng.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
        pop.N, &mu[0], bound_mmodels, bound_recmodels, interlocus_rec,
        std::bind(multilocus_additive(), std::placeholders::_1,
                  std::placeholders::_2, std::placeholders::_3),
        pop.neutral, pop.selected);
    if (!std::isfinite(wbar))
        {
            throw std::runtime_error("wbar not finite");
        }
}

BOOST_AUTO_TEST_SUITE_END()
