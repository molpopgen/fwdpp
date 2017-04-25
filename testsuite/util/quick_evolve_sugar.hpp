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
#include <fwdpp/interlocus_recombination.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/experimental/sample_diploid_mloc.hpp>
#include <testsuite/util/migpop.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>

template <typename singlepop_object_t>
void
simulate_singlepop(singlepop_object_t &pop, const unsigned simlen = 10,
                   const unsigned popsize = 5000)
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
                pop.mcounts, pop.N, popsize, 0.005,
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
            pop.N = popsize;
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * pop.N);
        }
}

template <typename singlepop_object_t, typename rng_type>
unsigned
simulate_singlepop(singlepop_object_t &pop, const rng_type &rng,
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

// Fitness function
struct multilocus_additive
{
    using result_type = double;
    inline double
    operator()(
        const multiloc_popgenmut_fixture::poptype::dipvector_t::value_type
            &diploid,
        const multiloc_popgenmut_fixture::poptype::gcont_t &gametes,
        const multiloc_popgenmut_fixture::poptype::mcont_t &mutations) const
    {
        using dip_t = multiloc_popgenmut_fixture::poptype::dipvector_t::
            value_type::value_type;
        return std::max(
            0., 1. + std::accumulate(diploid.begin(), diploid.end(), 0.,
                                     [&gametes, &mutations](const double d,
                                                            const dip_t &dip) {
                                         return d + KTfwd::additive_diploid()(
                                                        gametes[dip.first],
                                                        gametes[dip.second],
                                                        mutations)
                                                - 1.;
                                     }));
    }
};

template <typename poptype, typename rng_type, typename mmodel_vec,
          typename recmodel_vec, typename fitness_fxn>
inline unsigned
simulate_mlocuspop(poptype &pop, const rng_type &rng,
                   const mmodel_vec &mutmodels, const recmodel_vec &recmodels,
                   const fitness_fxn &fitness, const std::vector<double> &mu,
                   const std::vector<double> &rbw, unsigned &generation,
                   const unsigned simlen = 10)
/*!
  \brief Quick function for evolving a multilocus deme simulation
  \ingroup testing
  \note this version CAN be used on the same population object
 */
{
    unsigned g = generation;
    auto interlocus_rec = KTfwd::make_binomial_interlocus_rec(
        rng.get(), rbw.data(), rbw.size());
    for (; generation < g + simlen; ++generation)
        {
            double wbar = KTfwd::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, 1000, &mu[0], mutmodels, recmodels,
                interlocus_rec,
                std::bind(multilocus_additive(), std::placeholders::_1,
                          std::placeholders::_2, std::placeholders::_3),
                pop.neutral, pop.selected);
            assert(check_sum(pop.gametes, 8000));
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2000);
        }
    return g + simlen;
}

template <typename poptype, typename rng_type, typename mmodel_vec,
          typename recmodel_vec, typename fitness_fxn>
inline unsigned
simulate_mlocuspop_experimental(
    poptype &pop, const rng_type &rng, const mmodel_vec &mutmodels,
    const recmodel_vec &recmodels, const fitness_fxn &fitness,
    const std::vector<double> &mu, const std::vector<double> &rbw,
    unsigned &generation, const unsigned simlen = 10)
/*!
  \brief Quick function for evolving a multilocus deme simulation using
  KTfwd::experimental
  \ingroup testing
  \note this version CAN be used on the same population object
 */
{
    unsigned g = generation;
    auto interlocus_rec = KTfwd::make_binomial_interlocus_rec(
        rng.get(), rbw.data(), rbw.size());
    for (; generation < g + simlen; ++generation)
        {
            double wbar = KTfwd::experimental::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, 1000, &mu[0], mutmodels, recmodels,
                interlocus_rec,
                std::bind(multilocus_additive(), std::placeholders::_1,
                          std::placeholders::_2, std::placeholders::_3),
                pop.neutral, pop.selected);
            assert(check_sum(pop.gametes, 8000));
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2000);
        }
    return g + simlen;
}

template <typename metapop_object>
void
simulate_metapop(metapop_object &pop, const unsigned simlen = 10)
{
    // Evolve for 10 generations
    std::vector<std::function<double(
        const typename metapop_object::diploid_t &,
        const typename metapop_object::gcont_t &,
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
