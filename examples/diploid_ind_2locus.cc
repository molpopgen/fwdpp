/*! \include diploid_ind_2locus.cc
  Simple example of a two-locus simulation using the multilocus API in fwdpp.
*/

#include <iostream>
#include <fwdpp/diploid.hh>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include <vector>
#include <list>
#include <sstream>
// Use mutation model from sugar layer
#include <fwdpp/sugar/infsites.hpp>

using mtype = KTfwd::popgenmut;
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
    if (argc != 9)
        {
            std::cerr << "Incorrect number of arguments.\n"
                      << "Usage:\n"
                      << argv[0] << " N theta rho rbw ngens n nreps seed\n"
                      << "Where:\n"
                      << "N = population size (number of diploids)\n"
                      << "theta = 4Nu, the scaled neutral mutation rate\n"
                      << "rho = 4Nr, the scale recombination rate\n"
                      << "rbw = the probability that the two loci cross over, "
                         "per generation\n"
                      << "ngens = the number of generations to simulate\n"
                      << "n = the sample size to pull from the population at "
                         "the end of each simulated replicate\n"
                      << "nreps = the number of replicates to simulated\n"
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
    const unsigned ngens
        = atoi(argv[argument++]); // Number of generations to simulate
    const unsigned samplesize1
        = atoi(argv[argument++]); // Sample size to draw from the population
    int nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned seed = atoi(argv[argument++]); // Random number seed

    const std::vector<double> mu(
        2, theta / double(4 * N)); // per-gamete mutation rate per locus

    /*
      littler r is the recombination rate per region per generation.

      For individual simulation (UNLIKE GAMETE-BASED SIMS!!!),
      r = rho/(4N)
    */
    const double littler = rho / double(4 * N);

    // Initiate random number generation system
    GSLrng r(seed);

    unsigned twoN = 2 * N;

    std::function<double(void)> recmap = std::bind(gsl_rng_uniform, r.get()),
                                recmap2
                                = std::bind(gsl_ran_flat, r.get(), 1., 2.);

    while (nreps--)
        {
            multiloc_t pop(N, 2);
            assert(check_sum(pop.gametes, 2 * twoN));
            assert(popdata_sane_multilocus(pop.diploids, pop.gametes,
                                           pop.mutations, pop.mcounts));
            pop.mutations.reserve(size_t(
                2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
            unsigned generation = 0;
            double wbar;

            std::vector<std::function<std::vector<double>(
                const multiloc_t::gamete_t &, const multiloc_t::gamete_t &,
                const multiloc_t::mcont_t &)>>
                recpols{
                    std::bind(KTfwd::poisson_xover(), r.get(), littler, 0., 1.,
                              std::placeholders::_1, std::placeholders::_2,
                              std::placeholders::_3),
                    std::bind(KTfwd::poisson_xover(), r.get(), littler, 1., 2.,
                              std::placeholders::_1, std::placeholders::_2,
                              std::placeholders::_3)
                };

            std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                                  multiloc_t::mcont_t &)>>
                mmodels{
                    // Locus 0: positions Uniform [0,1)
                    std::bind(KTfwd::infsites(), std::placeholders::_1,
                              std::placeholders::_2, r.get(),
                              std::ref(pop.mut_lookup), &generation, mu[0], 0.,
                              [&r]() { return gsl_rng_uniform(r.get()); },
                              []() { return 0.; }, []() { return 0.; }),
                    // Locus 1: positions Uniform [1,2)
                    std::bind(KTfwd::infsites(), std::placeholders::_1,
                              std::placeholders::_2, r.get(),
                              std::ref(pop.mut_lookup), &generation, mu[1], 0.,
                              [&r]() { return gsl_ran_flat(r.get(), 1., 2.); },
                              []() { return 0.; }, []() { return 0.; })
                };

            std::vector<std::function<unsigned(void)>> interlocus_rec{
                std::bind(gsl_ran_binomial, r.get(), rbw, 1)};

            for (generation = 0; generation < ngens; ++generation)
                {
                    // Iterate the population through 1 generation
                    KTfwd::sample_diploid(
                        r.get(), pop.gametes, pop.diploids, pop.mutations,
                        pop.mcounts, N, &mu[0], mmodels, recpols,
                        interlocus_rec,
                        std::bind(no_selection_multi(), std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3),
                        pop.neutral, pop.selected);
                    assert(check_sum(pop.gametes, 2 * twoN));
                    KTfwd::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, 2 * N);
                    assert(popdata_sane_multilocus(pop.diploids, pop.gametes,
                                                   pop.mutations,
                                                   pop.mcounts));
                }
            // For giggles, make sure that the pop. is copy-constructible...
            multiloc_t pop2(pop);
            // Take a sample and print it to screen.
            auto x = KTfwd::ms_sample(r.get(), pop.mutations, pop.gametes,
                                      pop.diploids, samplesize1, true);
#ifdef HAVE_LIBSEQUENCE
            Sequence::SimData l1(x[0].begin(), x[0].end()),
                l2(x[1].begin(), x[1].end());
            std::cout << l1 << '\n' << l2 << '\n';
#endif
        }
    return 0;
}
