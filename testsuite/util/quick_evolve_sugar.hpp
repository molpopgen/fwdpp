/*!
 * \file quick_evolve_sugar.hpp
 * \brief Helper functions for unit/integration testing
 * \ingroup unit
 */
#ifndef FWDPP_TESTSUITE_UTIL_QUICK_EVOLVE_SUGAR_HPP
#define FWDPP_TESTSUITE_UTIL_QUICK_EVOLVE_SUGAR_HPP

#include <fwdpp/recombination.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/sample_diploid.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <testsuite/util/migpop.hpp>

template <typename singlepop_object_t>
void
simulate_singlepop(singlepop_object_t &pop, const unsigned simlen = 10)
/*!
  \brief Quick function for evolving a single-deme simulation
  \ingroup testing
  \note Do NOT call this function repeatedly on the same population.
 */
{
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);

    for (unsigned generation = 0; generation < simlen; ++generation)
        {
            double wbar = KTfwd::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, 1000, 0.005,
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, rng.get(),
                          std::ref(pop.mut_lookup), generation, 0.0025, 0.0025,
                          [&rng]() { return gsl_rng_uniform(rng.get()); },
                          []() { return -0.01; }, []() { return 1.; }),
                std::bind(KTfwd::poisson_xover(), rng.get(), 0.005, 0., 1.,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                std::bind(KTfwd::multiplicative_diploid(),
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, 2),
                pop.neutral, pop.selected);
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * pop.N);
        }
}

template <typename singlepop_object_t, typename rng_type>
unsigned
simulate_singlepop2(singlepop_object_t &pop, const rng_type &rng,
                    const unsigned generation, const unsigned simlen)
/*!
  \brief Quick function for evolving a single-deme simulation
  \ingroup testing
  \note this version CAN be used on the same population object
 */
{
    unsigned g = generation;
    for (; g < generation + simlen; ++g)
        {
            double wbar = KTfwd::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, 1000, 0.005,
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, rng.get(),
                          std::ref(pop.mut_lookup), g, 0.0025, 0.0025,
                          [&rng]() { return gsl_rng_uniform(rng.get()); },
                          []() { return -0.01; }, []() { return 1.; }),
                std::bind(KTfwd::poisson_xover(), rng.get(), 0.005, 0., 1.,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                std::bind(KTfwd::multiplicative_diploid(),
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, 2),
                pop.neutral, pop.selected);
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, g, 2 * pop.N);
        }
    return g + simlen;
}
template <typename metapop_object>
void
simulate_metapop(metapop_object &pop, const unsigned simlen = 10)
{
    // Evolve for 10 generations
    std::vector<std::function<double(
        const typename metapop_object::gamete_t &,
        const typename metapop_object::gamete_t &,
        const typename metapop_object::mcont_t &)>>
        fitness_funcs(2,
                      std::bind(KTfwd::multiplicative_diploid(),
                                std::placeholders::_1, std::placeholders::_2,
                                std::placeholders::_3, 2.));
    KTfwd::GSLrng_t<KTfwd::GSL_RNG_TAUS2> rng(0u);
    for (unsigned generation = 0; generation < simlen; ++generation)
        {
            std::vector<double> wbar = KTfwd::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, &pop.Ns[0], 0.005,
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, rng.get(),
                          std::ref(pop.mut_lookup), generation, 0.005, 0.,
                          [&rng]() { return gsl_rng_uniform(rng.get()); },
                          []() { return 0.; }, []() { return 0.; }),
                std::bind(KTfwd::poisson_xover(), rng.get(), 0.005, 0., 1.,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                fitness_funcs,
                std::bind(migpop, std::placeholders::_1, rng.get(), 0.001),
                pop.neutral, pop.selected);
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 4000);
        }
}

#endif
