/*! \include K_linked_regions_extensions.cc
  Simple example building up a K-locus simulation using the
  extensions sub-library
*/

#include <iostream>
#include <fwdpp/diploid.hh>
#include <vector>
#include <list>
#include <sstream>
// Use mutation model from sugar layer
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>

using mtype = KTfwd::popgenmut;
#define SINGLEPOP_SIM
#include <common_ind.hpp>

// This is our fitness model
struct additive_over_loci
{
    template <typename diploid_t, typename gcont_t, typename mcont_t>
    inline double
    operator()(const diploid_t &dip, const gcont_t &gametes,
               const mcont_t &mutations, const unsigned K) const noexcept
    /*
      The fitness model is additive over loci, multiplicative within loci
     */
    {
        double rv = 0.0;
        auto &g1 = gametes[dip.first];
        auto &g2 = gametes[dip.second];
        for (unsigned i = 0; i < K; ++i)
            {
                /*
                  Find the iterators corrsponding to locus 1 in g1 and g2.
                  Because gamete store mutation keys sorted by mutation
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

                rv += KTfwd::multiplicative_diploid()(start1, stop1, start2,
                                                      stop2, mutations, 2.0)
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
    const double rho = atof(argv[argument++]);   // 4*n*recombination rate.
                                                 // Note: recombination rate is
                                                 // per REGION, not SITE!!
    const double rbw = atof(argv[argument++]);   // rec rate b/w loci.
    const unsigned K = atoi(argv[argument++]); // Number of loci in simulation
    const unsigned ngens = atoi(argv[argument++]); //# generations to simulatae
    const unsigned seed
        = atoi(argv[argument++]); // Number of loci in simulation

    const double mutrate_region = theta / double(4 * N * K);
    const double mutrate_del_region = 0.1 * mutrate_region;
    const double recrate_region = rho / double(4 * N * K);

    std::vector<double> RBW(K, rbw);
    // Initiate random number generation system
    GSLrng r(seed);

    unsigned twoN = 2 * N;

    singlepop_t pop(N);
    pop.mutations.reserve(
        size_t(2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
    unsigned generation = 0;
    double wbar;

    // Set up mutation models
    std::vector<double> locus_starts(K);
    std::vector<double> locus_ends(K);
    for (unsigned i = 0; i < K; ++i)
        {
            locus_starts[i] = i;
            locus_ends[i] = i + 1;
        }
    std::vector<double> locus_weights(K, 1.0); // all loci will have equal
                                               // weight for both mutation and
                                               // recombination
    std::vector<KTfwd::extensions::shmodel> shmodels(
        K, KTfwd::extensions::shmodel(KTfwd::extensions::constant(-0.1),
                                      KTfwd::extensions::constant(1.0)));
    KTfwd::extensions::discrete_mut_model mmodels(
        locus_starts, locus_ends, locus_weights, locus_starts, locus_ends,
        locus_weights, shmodels);

    // Set up recombination rates
    locus_starts.clear();
    locus_ends.clear();
    locus_weights.clear();
    for (unsigned i = 0; i < K; ++i)
        {
            // start, end, and weight associated w/crossover w/in region i
            locus_starts.push_back(i);
            locus_ends.push_back(i + 1);
            locus_weights.push_back(recrate_region);

            // values associated w/crossover b/w region i and i+1
            locus_starts.push_back(i + 1);
            locus_ends.push_back(i + 1);
            locus_weights.push_back(rbw);
        }
    KTfwd::extensions::discrete_rec_model recmap(locus_starts, locus_ends,
                                                 locus_weights);

    // Now, synthesize a function that binds to the operator() and recmap

    const auto recpolicy = KTfwd::extensions::bind_drm(
        recmap, pop.gametes, pop.mutations, r.get(),
        double(K) * recrate_region + double(K - 1) * rbw);
    const auto bound_mmodels = KTfwd::extensions::bind_dmm(
        mmodels, pop.mutations, pop.mut_lookup, r.get(),
        double(K * mutrate_region), double(K * mutrate_del_region),
        &generation);

    for (generation = 0; generation < ngens; ++generation)
        {
            // Iterate the population through 1 generation
            KTfwd::sample_diploid(
                r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                N, double(K) * (mutrate_region + mutrate_del_region),
                // This is the synthesized function bound to operator() of
                // mmodels:
                bound_mmodels, recpolicy,
                std::bind(additive_over_loci(), std::placeholders::_1,
                          std::placeholders::_2, std::placeholders::_3, K),
                pop.neutral, pop.selected);
            assert(check_sum(pop.gametes, K * twoN));
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
            assert(popdata_sane(pop.diploids, pop.gametes, pop.mutations,
                                pop.mcounts));
        }
}
