/*! \include K_linked_regions_generalized_rec.cc
  Simple example building up variation in recombination
  rate using fwdpp::general_rec_variation
*/

#include <iostream>
#include <fwdpp/diploid.hh>
#include <fwdpp/recbinder.hpp>
#ifdef HAVE_LIBSEQUENCE
#endif
#include <vector>
#include <list>
#include <sstream>
// Use mutation model from sugar layer
#include <fwdpp/popgenmut.hpp>
using mtype = fwdpp::popgenmut;
#define DIPLOID_POPULATION_SIM
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

    const double mu = theta / double(4 * N);

    /*
      littler r is the recombination rate per region per generation.
    */
    const double littler = rho / double(4 * N * K);

    // Initiate random number generation system
    GSLrng r(seed);

    diploid_population pop(N);
    pop.mutations.reserve(
        size_t(2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
    unsigned generation = 0;
    double wbar;

    fwdpp::general_rec_variation recvar;
    for (unsigned i = 0; i < K; ++i)
        {
            recvar.recmap.push_back(
                fwdpp::poisson_interval(littler, i, i + 1.0));
            if (i)
                {
                    recvar.recmap.push_back(
                        fwdpp::crossover_point(rbw, i + 1.0, false));
                }
        }
    auto recmap = fwdpp::recbinder(recvar, r.get());

    const auto mmodel
        = [&pop, &r, &generation, K](fwdpp::flagged_mutation_queue &recbin,
                                     diploid_population::mcont_t &mutations) {
              return fwdpp::infsites_popgenmut(
                  recbin, mutations, r.get(), pop.mut_lookup, generation, 0.0,
                  [&r, K]() { return gsl_ran_flat(r.get(), 0, K); },
                  []() { return 0.0; }, []() { return 0.0; });
          };

    for (generation = 0; generation < ngens; ++generation)
        {
            // Iterate the population through 1 generation
            fwdpp::sample_diploid(
                r.get(), pop.haploid_genomes, pop.diploids, pop.mutations,
                pop.mcounts, N, mu, mmodel, recmap,
                fwdpp::multiplicative_diploid(fwdpp::fitness(2.)), pop.neutral,
                pop.selected);
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
            fwdpp::debug::validate_sum_haploid_genome_counts(
                pop.haploid_genomes, 2 * N);
            fwdpp::debug::validate_pop_data(pop);
        }
#ifdef HAVE_LIBSEQUENCE
#endif
}
