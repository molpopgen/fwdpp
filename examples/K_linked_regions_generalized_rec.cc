/*! \include K_linked_regions_generalized_rec.cc
  Simple example building up variation in recombination
  rate using KTfwd::general_rec_variation
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
#include <fwdpp/sugar/sampling.hpp>
using mtype = KTfwd::popgenmut;
#define SINGLEPOP_SIM
#include <common_ind.hpp>

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

	const double mu = theta/double(4*N);

    /*
      littler r is the recombination rate per region per generation.
    */
    const double littler = rho / double(4 * N * K);

    // Initiate random number generation system
    GSLrng r(seed);

    unsigned twoN = 2 * N;

    singlepop_t pop(N);
    pop.mutations.reserve(
        size_t(2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
    unsigned generation = 0;
    double wbar;

    KTfwd::general_rec_variation recvar;
    for (unsigned i = 0; i < K; ++i)
        {
            recvar.recmap.push_back(
                KTfwd::poisson_interval(r.get(), littler, i, i + 1.0));
            if (i)
                {
                    recvar.recmap.push_back(
                        KTfwd::crossover_point(r.get(), rbw, i + 1.0, false));
                }
        }

    for (generation = 0; generation < ngens; ++generation)
        {
            // Iterate the population through 1 generation
            KTfwd::sample_diploid(
                r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                N,mu, std::bind(KTfwd::infsites(), std::placeholders::_1,
                             std::placeholders::_2, r.get(),
                             std::ref(pop.mut_lookup), generation, mu, 0.,
                             [&r,K]() { return gsl_ran_flat(r.get(),0.,double(K)); },
                             []() { return 0.; }, []() { return 0.; }),
                recvar, std::bind(KTfwd::multiplicative_diploid(),
                                  std::placeholders::_1, std::placeholders::_2,
                                  std::placeholders::_3, 2.),
                pop.neutral, pop.selected);
            assert(check_sum(pop.gametes, K * twoN));
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
        }
    auto x = KTfwd::ms_sample(r.get(), pop.mutations, pop.gametes,
                              pop.diploids, 10, true);
#ifdef HAVE_LIBSEQUENCE
    Sequence::SimData a(x.begin(), x.end());
    std::cout << a << '\n';
#endif
}
