/*!
 * \file quick_evolve_sugar.hpp
 * \brief Helper functions for unit/integration testing
 * \ingroup unit
 */
#ifndef FWDPP_TESTSUITE_UTIL_QUICK_EVOLVE_SUGAR_HPP
#define FWDPP_TESTSUITE_UTIL_QUICK_EVOLVE_SUGAR_HPP

#include <exception>
#include <cmath>
#include <fwdpp/debug.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/sample_diploid.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/interlocus_recombination.hpp>
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/GSLrng_t.hpp>
#include <testsuite/fixtures/sugar_fixtures.hpp>

template <typename slocuspop_object_t>
void
simulate_slocuspop(slocuspop_object_t &pop, const unsigned simlen = 10,
                   const unsigned popsize = 5000)
/*!
  \brief Quick function for evolving a single-deme simulation
  \ingroup testing
  \note Do NOT call this function repeatedly on the same population.
 */
{
    fwdpp::GSLrng_t<fwdpp::GSL_RNG_TAUS2> rng(0u);
    unsigned generation = 0;
    const auto mmodel = [&pop, &rng, &generation](
                            fwdpp::flagged_mutation_queue &recbin,
                            typename slocuspop_object_t::mcont_t &mutations) {
        return fwdpp::infsites_popgenmut(
            recbin, mutations, rng.get(), pop.mut_lookup, generation, 0.5,
            [&rng]() { return gsl_rng_uniform(rng.get()); },
            []() { return -0.01; }, []() { return 1.; });
    };
    for (; generation < simlen; ++generation)
        {
            double wbar = fwdpp::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, pop.N, popsize, 0.005, mmodel,
                fwdpp::recbinder(fwdpp::poisson_xover(0.005, 0., 1.),
                                 rng.get()),
                fwdpp::multiplicative_diploid(fwdpp::fitness(2.)), pop.neutral,
                pop.selected);
            if (!std::isfinite(wbar))
                {
                    throw std::runtime_error("fitness not finite");
                }
            pop.N = popsize;
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * pop.N);
        }
}

template <typename slocuspop_object_t, typename rng_type>
unsigned
simulate_slocuspop(slocuspop_object_t &pop, const rng_type &rng,
                   const unsigned generation, const unsigned simlen)
/*!
  \brief Quick function for evolving a single-deme simulation
  \ingroup testing
  \note this version CAN be used on the same population object
 */
{
    unsigned g = generation;

    const auto mmodel = [&pop, &rng, &generation](
                            fwdpp::flagged_mutation_queue &recbin,
                            typename slocuspop_object_t::mcont_t &mutations) {
        return fwdpp::infsites_popgenmut(
            recbin, mutations, rng.get(), pop.mut_lookup, generation, 0.5,
            [&rng]() { return gsl_rng_uniform(rng.get()); },
            []() { return -0.01; }, []() { return 1.; });
    };
    for (; g < generation + simlen; ++g)
        {
            double wbar = fwdpp::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, 1000, 0.005, mmodel,
                fwdpp::recbinder(fwdpp::poisson_xover(0.005, 0., 1.),
                                 rng.get()),
                fwdpp::multiplicative_diploid(2.), pop.neutral, pop.selected);
            if (!std::isfinite(wbar))
                {
                    throw std::runtime_error("fitness not finite");
                }
            fwdpp::update_mutations(pop.mutations, pop.fixations,
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
        const mlocuspop_popgenmut_fixture::poptype::dipvector_t::value_type
            &diploid,
        const mlocuspop_popgenmut_fixture::poptype::gcont_t &gametes,
        const mlocuspop_popgenmut_fixture::poptype::mcont_t &mutations) const
    {
        using dip_t = mlocuspop_popgenmut_fixture::poptype::dipvector_t::
            value_type::value_type;
        return std::max(
            0., 1.
                    + std::accumulate(
                          diploid.begin(), diploid.end(), 0.,
                          [&gametes, &mutations](const double d,
                                                 const dip_t &dip) {
                              return d
                                     + fwdpp::additive_diploid(fwdpp::fitness(
                                           1.))(gametes[dip.first],
                                                gametes[dip.second], mutations)
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
    auto interlocus_rec = fwdpp::make_binomial_interlocus_rec(
        rng.get(), rbw.data(), rbw.size());
    for (; generation < g + simlen; ++generation)
        {
            double wbar = fwdpp::sample_diploid(
                rng.get(), pop.gametes, pop.diploids, pop.mutations,
                pop.mcounts, 1000, &mu[0], mutmodels, recmodels,
                interlocus_rec, fitness, pop.neutral, pop.selected);
            if (!std::isfinite(wbar))
                {
                    throw std::runtime_error("fitness not finite");
                }
            fwdpp::debug::validate_sum_gamete_counts(pop.gametes, 8000);
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2000);
        }
    return g + simlen;
}

#endif
