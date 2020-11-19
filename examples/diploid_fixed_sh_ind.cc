/*
  \include diploid_fixed_sh_ind.cc

  Simulate a recombining region with neutral mutations and non-neutral
  mutations with fixed 's' and 'h'.
 */

#include <fwdpp/recbinder.hpp>
#include <numeric>
#include <functional>
#include <cassert>
#include <iomanip>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/sample_diploid.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/debug.hpp>
#include <fwdpp/types/mutation.hpp>
#include <fwdpp/genetic_map/genetic_map.hpp>
#include <fwdpp/genetic_map/poisson_interval.hpp>
#include <fwdpp/algorithm/compact_mutations.hpp>
#define DIPLOID_POPULATION_SIM
// the type of mutation
using mtype = fwdpp::mutation;
#include <common_ind.hpp>
#include <gsl/gsl_randist.h>

int
main(int argc, char **argv)
{
    if (argc != 11)
        {
            std::cerr << "Too few arguments.\n"
                      << "Usage: diploid_fixed_sh_ind N theta_neutral "
                         "theta_deleterious rho s h ngens samplesize nreps "
                         "seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N = atoi(argv[argument++]);
    const double theta_neutral = atof(argv[argument++]);
    const double theta_del = atof(argv[argument++]);
    const double rho = atof(argv[argument++]);
    const double s = atof(argv[argument++]);
    const double h = atof(argv[argument++]);
    const unsigned ngens = atoi(argv[argument++]);
    const unsigned samplesize1 = atoi(argv[argument++]);
    int nreps = atoi(argv[argument++]);
    const unsigned seed = atoi(argv[argument++]);

    const double mu_neutral = theta_neutral / double(4 * N);
    const double mu_del = theta_del / double(4 * N);
    const double littler = rho / double(4 * N);

    std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    GSLrng r(seed);

    // recombination map is uniform[0,1)
    fwdpp::genetic_map gmap;
    gmap.add_callback(fwdpp::poisson_interval(0, 1, littler));
    const auto rec = fwdpp::recbinder(std::cref(gmap), r.get());
    const double pselected = mu_del / (mu_del + mu_neutral);

    while (nreps--)
        {
            diploid_population pop(N);
            pop.mutations.reserve(
                size_t(std::ceil(std::log(2 * N) * (theta_neutral + theta_del)
                                 + 0.667 * (theta_neutral + theta_del))));
            unsigned generation = 0;
            const auto mmodel = [&pop, &r, &generation, s, h, pselected](
                                    fwdpp::flagged_mutation_queue &recbin,
                                    diploid_population::mutation_container &mutations) {
                return fwdpp::infsites_mutation(
                    recbin, mutations, r.get(), pop.mut_lookup, generation, pselected,
                    [&r]() { return gsl_rng_uniform(r.get()); }, [s]() { return s; },
                    [h]() { return h; });
            };

            double wbar = 1;
            for (generation = 0; generation < ngens; ++generation)
                {
                    wbar = fwdpp::sample_diploid(
                        r.get(), pop.haploid_genomes, pop.diploids, pop.mutations,
                        pop.mcounts, N, mu_neutral + mu_del, mmodel,
                        // The function to generation recombination positions:
                        rec, fwdpp::multiplicative_diploid(fwdpp::fitness(1.)),
                        pop.neutral, pop.selected);
                    fwdpp::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, 2 * N);
                    if (generation && generation % 100 == 0.0)
                        {
                            fwdpp::compact_mutations(pop);
                        }
                    fwdpp::debug::validate_sum_haploid_genome_counts(pop.haploid_genomes,
                                                                     2 * N);
                }
            for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i])
                        {
                            std::cout << pop.mutations[i].s << ' ' << pop.mcounts[i]
                                      << '\n';
                        }
                }
            // Take a sample of size samplesize1.  Two data blocks are
            // returned, one for neutral mutations, and one for selected
            std::vector<std::size_t> random_dips;
            for (unsigned i = 0; i < samplesize1; ++i)
                {
                    auto x = static_cast<std::size_t>(gsl_ran_flat(r.get(), 0, N));
                    while (std::find(random_dips.begin(), random_dips.end(), x)
                           != random_dips.end())
                        {
                            x = static_cast<std::size_t>(gsl_ran_flat(r.get(), 0, N));
                        }
                }
            auto dm = fwdpp::sample_individuals(pop, random_dips, true, false, true);
        }
}
