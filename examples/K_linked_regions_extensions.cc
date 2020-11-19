/*! \include K_linked_regions_extensions.cc
  Simple example building up a K-locus simulation using the
  extensions sub-library
*/

#include <iostream>
#include <vector>
#include <list>
#include <sstream>
#include <fwdpp/debug.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/sample_diploid.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/util.hpp>
// Use mutation model from sugar layer
#include <fwdpp/types/mutation.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>

using mtype = fwdpp::mutation;
#define DIPLOID_POPULATION_SIM
#include <common_ind.hpp>

// This is our fitness model
struct additive_over_loci
{
    template <typename diploid_t, typename gcont_t, typename mcont_t>
    inline double
    operator()(const diploid_t &dip, const gcont_t &haploid_genomes,
               const mcont_t &mutations, const unsigned K) const noexcept
    /*
      The fitness model is additive over loci, multiplicative within loci
     */
    {
        double rv = 0.0;
        auto &g1 = haploid_genomes[dip.first];
        auto &g2 = haploid_genomes[dip.second];
        for (unsigned i = 0; i < K; ++i)
            {
                /*
                  Find the iterators corrsponding to locus 1 in g1 and g2.
                  Because haploid_genome store mutation keys sorted by mutation
                  position,
                  we can use efficient binary searches.

                  We want the positions of all mutations such that i <= p <
                  i+1,
                  which correspond to lower_bound and upper_bound,
                  respectively.
                */
                auto stop1 = std::upper_bound(
                    g1.smutations.begin(), g1.smutations.end(), double(i + 1),
                    [&mutations](const double val, const std::size_t i) {
                        return val < mutations[i].pos;
                    });
                auto start1 = std::lower_bound(
                    g1.smutations.begin(), stop1, double(i),
                    [&mutations](const std::size_t i, const double val) {
                        return mutations[i].pos < val;
                    });
                auto stop2 = std::upper_bound(
                    g2.smutations.begin(), g2.smutations.end(), double(i + 1),
                    [&mutations](const double val, const std::size_t i) {
                        return val < mutations[i].pos;
                    });
                auto start2 = std::lower_bound(
                    g2.smutations.begin(), stop2, double(i),
                    [&mutations](const std::size_t i, const double val) {
                        return mutations[i].pos < val;
                    });

                rv += fwdpp::multiplicative_diploid(fwdpp::fitness(2.0))(
                          start1, stop1, start2, stop2, mutations)
                      - 1.0;
            }
        return std::max(1.0 + rv, 0.0);
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
    const double rbw = atof(argv[argument++]);     // rec rate b/w loci.
    const unsigned K = atoi(argv[argument++]);     // Number of loci in simulation
    const unsigned ngens = atoi(argv[argument++]); //# generations to simulatae
    const unsigned seed = atoi(argv[argument++]);  // Number of loci in simulation

    const double mutrate_region = theta / double(4 * N * K);
    const double mutrate_del_region = 0. * mutrate_region;
    const double recrate_region = rho / double(4 * N * K);

    std::vector<double> RBW(K, rbw);
    // Initiate random number generation system
    GSLrng r(seed);

    diploid_population pop(N);
    pop.mutations.reserve(
        size_t(2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
    unsigned generation = 0;
    double wbar;

    // Set up mutation models
    std::vector<double> locus_starts(K);
    std::vector<double> locus_ends(K);
    std::vector<double> locus_weights, locus_rec_weights;
    std::vector<fwdpp::traits::mutation_model<diploid_population::mutation_container>>
        functions;
    std::vector<fwdpp::extensions::discrete_rec_model::function_type> rec_functions;

    const double pselected = mutrate_del_region / (mutrate_del_region + mutrate_region);
    for (unsigned i = 0; i < K; ++i)
        {
            locus_weights.push_back(1.0);
            functions.push_back([&pop, &r, &generation, pselected,
                                 i](fwdpp::flagged_mutation_queue &recbin,
                                    diploid_population::mutation_container &mutations) {
                return fwdpp::infsites_mutation(
                    recbin, mutations, r.get(), pop.mut_lookup, generation, pselected,
                    [&r, i]() { return gsl_ran_flat(r.get(), i, i + 1); },
                    []() { return 0.0; }, []() { return 0.0; });
            });
            rec_functions.push_back([i, &r](std::vector<double> &breakpoints) {
                breakpoints.push_back(gsl_ran_flat(r.get(), i, i + 1));
            });
            locus_rec_weights.push_back(recrate_region);
            if (i < K - 1)
                {
                    rec_functions.push_back([i](std::vector<double> &breakpoints) {
                        breakpoints.push_back(i + 1);
                    });
                    locus_rec_weights.push_back(rbw);
                }
        }

    fwdpp::extensions::discrete_mut_model<diploid_population::mutation_container>
        mmodels(std::move(functions), std::move(locus_weights));

    const auto bound_mmodels = fwdpp::extensions::bind_dmm(r.get(), mmodels);

    const double ttl_recrate = double(K) * recrate_region + double(K - 1) * rbw;

    fwdpp::extensions::discrete_rec_model recmap(ttl_recrate, std::move(rec_functions),
                                                 std::move(locus_rec_weights));

    auto bound_recmap = [&recmap, &r]() { return recmap(r.get()); };

    for (generation = 0; generation < ngens; ++generation)
        {
            // Iterate the population through 1 generation
            fwdpp::sample_diploid(
                r.get(), pop.haploid_genomes, pop.diploids, pop.mutations, pop.mcounts,
                N, double(K) * (mutrate_region + mutrate_del_region),
                // This is the synthesized function bound to operator() of
                // mmodels:
                bound_mmodels, bound_recmap,
                std::bind(additive_over_loci(), std::placeholders::_1,
                          std::placeholders::_2, std::placeholders::_3, K),
                pop.neutral, pop.selected);
            fwdpp::update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                                    pop.mut_lookup, pop.mcounts, generation, 2 * N);
            fwdpp::debug::validate_sum_haploid_genome_counts(pop.haploid_genomes, 2 * N);
            fwdpp::debug::validate_pop_data(pop);
        }
    for (unsigned i = 0; i < pop.mcounts.size(); ++i)
        {
            if (pop.mcounts[i])
                {
                    std::cout << pop.mutations[i].pos << ' ' << pop.mutations[i].s << ' '
                              << pop.mutations[i].g << ' ' << pop.mcounts[i] << '\n';
                }
        }
}
