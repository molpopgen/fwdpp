/*! \include K_linked_regions_multilocus.cc
  Simple example building up a K-locus simulation using the
  mulitlocus API
*/

#include <iostream>
#include <fwdpp/diploid.hh>
#include <fwdpp/recbinder.hpp>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include <vector>
#include <list>
#include <sstream>
// Use mutation model from sugar layer
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/sampling.hpp>
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

    multiloc_t pop(N, K);
    pop.mutations.reserve(
        size_t(2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
    unsigned generation = 0;
    double wbar;

    std::vector<std::function<std::vector<double>()>> recpols;
    std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                          multiloc_t::mcont_t &)>>
        mmodels;
    for (unsigned i = 0; i < K; ++i)
        {
            recpols.emplace_back(fwdpp::recbinder(
                fwdpp::poisson_xover(littler, i, i + 1), r.get()));
            mmodels.push_back(
                [&pop, &r, &generation](std::queue<std::size_t> &recbin,
                                        multiloc_t::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, r.get(), pop.mut_lookup, generation,
                        0.0, [&r]() { return gsl_rng_uniform(r.get()); },
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
            assert(check_sum(pop.gametes, K * twoN));
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
            assert(popdata_sane_multilocus(pop.diploids, pop.gametes,
                                           pop.mutations, pop.mcounts));
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
    auto x = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                              pop.diploids, 10, true);
#ifdef HAVE_LIBSEQUENCE
    for (auto &i : x)
        {
            Sequence::SimData a(i.begin(), i.end());
            std::cout << a << '\n';
        }
#endif
}
