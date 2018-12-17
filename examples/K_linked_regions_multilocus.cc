/*! \include K_linked_regions_multilocus.cc
  Simple example building up a K-locus simulation using the
  mulitlocus API
*/

#include <iostream>
#include <fwdpp/diploid.hh>
#include <fwdpp/recbinder.hpp>
#ifdef HAVE_LIBSEQUENCE
#endif
#include <vector>
#include <list>
#include <sstream>
#include <fwdpp/debug.hpp>
// Use mutation model from sugar layer
#include <fwdpp/popgenmut.hpp>
#include <fwdpp/sampling_functions.hpp>
using mtype = fwdpp::popgenmut;
#define MULTILOCUS_SIM
#include <common_ind.hpp>

// Fitness function
struct no_selection_multi
{
    typedef double result_type;
    inline double
    operator()(const multiloc_t::dipvector_t::value_type &,
               const multiloc_t::gcont_t &, const multiloc_t::mcont_t &) const
    {
        return 1.;
    }
};

int
main(int argc, char **argv)
{
    int argument = 1;
    if (argc != 8)
        {
            std::cerr << "Incorrect number of arguments.\n"
                      << "Usage:\n"
                      << argv[0] << " N theta rho rbw K ngens seed\n"
                      << "Where:\n"
                      << "N = population size (number of diploids)\n"
                      << "theta = 4Nu, the scaled neutral mutation "
                         "rate,summed over all loci\n"
                      << "rho = 4Nr, the scale recombination rate, summed "
                         "over all loci\n"
                      << "rbw = the recombination rate between locus i and "
                         "i+1, constant for all i.\n"
                      << "K = the number of loci in the simulation\n"
                      << "ngens = length of simulation in generatios\n"
                      << "seed = seed value for random number generations\n";
            std::exit(0);
        }
    const unsigned N = atoi(argv[argument++]);   // Number of diploids
    const double theta = atof(argv[argument++]); // 4*n*mutation rate.  Note:
    // mutation rate is per
    // REGION, not SITE!!
    const double rho = atof(argv[argument++]); // 4*n*recombination rate.
    // Note: recombination rate is
    // per REGION, not SITE!!
    const double rbw = atof(argv[argument++]); // rec rate b/w loci.
    const unsigned K = atoi(argv[argument++]); // Number of loci in simulation
    const unsigned ngens = atoi(argv[argument++]); //# generations to simulatae
    const unsigned seed
        = atoi(argv[argument++]); // Number of loci in simulation

    const std::vector<double> mu(
        K, theta / double(4 * N * K)); // per-gamete mutation rate per locus

    /*
      littler r is the recombination rate per region per generation.
    */
    const double littler = rho / double(4 * N * K);

    // Initiate random number generation system
    GSLrng r(seed);

    std::vector<std::pair<double, double>> boundaries;
    for (int i = 0; i < K; ++i)
        {
            boundaries.emplace_back(i, i + 1);
        }
    multiloc_t pop(N, boundaries);
    pop.mutations.reserve(
        size_t(2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
    unsigned generation = 0;

    std::vector<std::function<std::vector<double>()>> recpols;
    std::vector<std::function<std::size_t(fwdpp::flagged_mutation_queue &,
                                          multiloc_t::mcont_t &)>>
        mmodels;
    for (unsigned i = 0; i < K; ++i)
        {
            pop.locus_boundaries.emplace_back(i, i + 1);
            recpols.emplace_back(fwdpp::recbinder(
                fwdpp::poisson_xover(littler, i, i + 1), r.get()));
            mmodels.push_back([&pop, &r, &generation,
                               i](fwdpp::flagged_mutation_queue &recbin,
                                  multiloc_t::mcont_t &mutations) {
                return fwdpp::infsites_popgenmut(
                    recbin, mutations, r.get(), pop.mut_lookup, generation,
                    0.0, [&r, i]() { return gsl_ran_flat(r.get(), i, i + 1); },
                    []() { return 0.0; }, []() { return 0.0; });
            });
        }
    std::vector<std::function<unsigned(void)>> interlocus_rec(
        K - 1, std::bind(gsl_ran_binomial, r.get(), rbw, 1));
    for (generation = 0; generation < ngens; ++generation)
        {
            // Iterate the population through 1 generation
            fwdpp::sample_diploid(
                r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                N, mu.data(), mmodels, recpols, interlocus_rec,
                no_selection_multi(), pop.neutral, pop.selected);
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
            fwdpp::debug::validate_sum_gamete_counts(
                pop.gametes, K * 2 * pop.diploids.size());
            fwdpp::debug::validate_pop_data(pop);
#ifndef NDEBUG
            /*
            Useful block for long-run testing.

            Make sure that all mutation positions are in the "right" locus.
                 */
            for (const auto &dip : pop.diploids)
                {
                    for (unsigned locus = 0; locus < dip.size(); ++locus)
                        {
                            assert(!std::any_of(
                                pop.gametes[dip[locus].first]
                                    .mutations.begin(),
                                pop.gametes[dip[locus].first].mutations.end(),
                                [&pop, locus](const std::size_t mi) {
                                    return !(pop.mutations[mi].pos
                                                 >= double(locus)
                                             && pop.mutations[mi].pos
                                                    < double(locus) + 1.0);
                                }));
                        }
                }
#endif
        }
    // Take a sample of size samplesize1 from the population
    std::vector<std::size_t> random_dips;
    for (unsigned i = 0; i < N / 100; ++i)
        {
            auto x = static_cast<std::size_t>(gsl_ran_flat(r.get(), 0, N));
            while (std::find(random_dips.begin(), random_dips.end(), x)
                   != random_dips.end())
                {
                    x = static_cast<std::size_t>(gsl_ran_flat(r.get(), 0, N));
                }
            random_dips.push_back(x);
        }
    auto s = fwdpp::sample_individuals_by_window(
        pop, random_dips, pop.locus_boundaries, true, true, true);
#ifdef HAVE_LIBSEQUENCE
#endif
}
