/*! \include K_linked_regions_multilocus.cc
  Simple example building up a K-locus simulation using the
  mulitlocus API
*/

#include <iostream>
#include <fwdpp/diploid.hh>
#include <Sequence/SimData.hpp>
#include <vector>
#include <list>
#include <sstream>
// Use mutation model from sugar layer
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/sampling.hpp>
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

    std::vector<double> RBW(K, rbw);
    // Initiate random number generation system
    GSLrng r(seed);

    unsigned twoN = 2 * N;

    multiloc_t pop(N, K);
    pop.mutations.reserve(
        size_t(2 * std::ceil(std::log(2 * N) * (theta) + 0.667 * (theta))));
    unsigned generation = 0;
    double wbar;

    std::vector<std::function<std::vector<double>(
        const multiloc_t::gamete_t &, const multiloc_t::gamete_t &,
        const multiloc_t::mcont_t &)>>
        recpols;
    std::vector<std::function<std::size_t(std::queue<std::size_t> &,
                                          multiloc_t::mcont_t &)>>
        mmodels;
    for (unsigned i = 0; i < K; ++i)
        {
            recpols.push_back(
                std::bind(KTfwd::poisson_xover(), r.get(), littler, double(i),
                          double(i) + 1.0, std::placeholders::_1,
                          std::placeholders::_2, std::placeholders::_3));
            mmodels.push_back(std::bind(
                KTfwd::infsites(), std::placeholders::_1,
                std::placeholders::_2, r.get(), std::ref(pop.mut_lookup),
                &generation, mu[i], 0.,
                [&r, i]() {
                    return gsl_ran_flat(r.get(), double(i), double(i) + 1.0);
                },
                []() { return 0.; }, []() { return 0.; }));
        }

    for (generation = 0; generation < ngens; ++generation)
        {
            // Iterate the population through 1 generation
            KTfwd::sample_diploid(
                r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                N, mu.data(), mmodels, recpols, RBW.data(),
                [](const gsl_rng *__r, const double __d) -> unsigned {
                    return (gsl_rng_uniform(__r) <= __d) ? 1u : 0u;
                },
                std::bind(no_selection_multi(), std::placeholders::_1,
                          std::placeholders::_2, std::placeholders::_3),
                pop.neutral, pop.selected);
            assert(check_sum(pop.gametes, K * twoN));
            KTfwd::update_mutations(pop.mutations, pop.fixations,
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
    auto x = KTfwd::ms_sample(r.get(), pop.mutations, pop.gametes,
                              pop.diploids, 10, true);
    for (auto &f : pop.fixations)
        std::cout << f.pos << ' ';
    std::cout << '\n';
    for (auto &i : x)
        {
            Sequence::SimData a(i.begin(), i.end());
            std::cout << a << '\n';
        }
    std::vector<std::pair<double, double>> boundaries;
    for (unsigned i = 0; i < K; ++i)
        boundaries.emplace_back(i, i + 1);
    auto y = KTfwd::sample(r.get(), pop, 10, false, boundaries);
    for (auto &&i : y)
        {
            Sequence::SimData a(i.begin(), i.end());
            std::cout << a << '\n';
        }
}
