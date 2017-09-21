/*
  \include bneck_selection.cc

  Bottleneck + exponential recovery.  With selection.
 */

#include <numeric>
#include <cmath>
#include <functional>
#include <cassert>
#include <iomanip>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include <fwdpp/diploid.hh>
// Pull mutation model from fwdpp's "sugar" layer  (@ref md_md_sugar)
#include <fwdpp/sugar/infsites.hpp>

// typedef mutation_with_age mtype;
using mtype = KTfwd::mutation;
#define SINGLEPOP_SIM
#include <common_ind.hpp>

using namespace KTfwd;

int
main(int argc, char **argv)
{
    if (argc != 14)
        {
            std::cerr << "Error, too few arguments.\n"
                      << "Usage: " << argv[0] << ' '
                      << "N theta_neutral theta_deleterious 4Nr s h ngens N2 "
                         "N3 ngens2 n nreps seed\n";
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
    const unsigned N2
        = atoi(argv[argument++]); // change N to N2 after ngens of evolution
    const unsigned N3 = atoi(
        argv[argument++]); // N2 will change to N2 during ngens2 of exp. growth
    const unsigned ngens2 = atoi(argv[argument++]);
    const unsigned samplesize1 = atoi(argv[argument++]);
    int nreps = atoi(argv[argument++]);
    const unsigned seed = atoi(argv[argument++]);

    const double mu_neutral = theta_neutral / double(4 * N);
    const double mu_del = theta_del / double(4 * N);
    const double littler = rho / double(4 * N);

    // Do some basic argument checking
    if (N2 > N)
        {
            std::cerr << "Error, N2 > N (" << N2 << " > " << N
                      << "), but it should be N2 <= N\n";
        }
    if (N2 > N3)
        {
            std::cerr << "Error, N3 > N2 (" << N3 << " > " << N2
                      << "), but it should be N3 > N2\n";
        }
    if (ngens2 == 0)
        {
            std::cerr << "Error, ngens2 equals zero.  Must be > 0\n";
        }
    std::copy(argv, argv + argc,
              std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    GSLrng r(seed);

    // recombination map is uniform[0,1)
    std::function<double(void)> recmap = std::bind(gsl_rng_uniform, r.get());

    std::vector<unsigned> N_over_time(ngens, N + 1);
    N_over_time.push_back(N2);
    // Figure out the growth rate, etc.
    double G = std::exp((std::log(double(N3)) - std::log(double(N2)))
                        / double(ngens2));
    for (unsigned i = 0; i < ngens2; ++i)
        {
            N_over_time.push_back(unsigned(round(N2 * std::pow(G, i + 1))));
        }

    while (nreps--)
        {
            // Initialize a population of N diploids via KTfwd::singlepop
            // (fwdpp/sugar/singlepop.hpp)
            singlepop_t pop(N);
            pop.mutations.reserve(
                std::ceil(std::log(2 * N) * (theta_neutral + theta_del)
                          + 0.667 * (theta_neutral + theta_del)));
            unsigned generation = 0;
            double wbar = 1;
            for (unsigned generation = 0; generation < N_over_time.size() - 1;
                 ++generation)
                {
                    assert(KTfwd::check_sum(pop.gametes, 2 * N));
                    wbar = KTfwd::sample_diploid(
                        r.get(), pop.gametes, pop.diploids, pop.mutations,
                        pop.mcounts, N_over_time[generation],
                        N_over_time[generation + 1], mu_neutral + mu_del,
                        std::bind(KTfwd::infsites(), std::placeholders::_1,
                                  std::placeholders::_2, r.get(),
                                  std::ref(pop.mut_lookup), mu_neutral, mu_del,
                                  [&r]() { return gsl_rng_uniform(r.get()); },
                                  [&s]() { return s; }, [&h]() { return h; }),
                        // The function to generation recombination positions:
                        std::bind(KTfwd::poisson_xover(), r.get(), littler, 0.,
                                  1., std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3),
                        std::bind(KTfwd::multiplicative_diploid(),
                                  std::placeholders::_1, std::placeholders::_2,
                                  std::placeholders::_3, 2.),
                        pop.neutral, pop.selected);
                    KTfwd::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, 2 * N);
                    assert(KTfwd::check_sum(pop.gametes, 2 * N));
                }

            // Take a sample of size samplesize1.  Two data blocks are
            // returned, one for neutral mutations, and one for selected
            std::pair<std::vector<std::pair<double, std::string>>,
                      std::vector<std::pair<double, std::string>>>
                sample
                = ms_sample_separate(r.get(), pop.mutations, pop.gametes,
                                     pop.diploids, samplesize1);

#ifdef HAVE_LIBSEQUENCE
            Sequence::SimData neutral_muts, selected_muts;
            neutral_muts.assign(sample.first.begin(), sample.first.end());
            selected_muts.assign(sample.second.begin(), sample.second.end());

            std::cout << neutral_muts << '\n' << selected_muts << '\n';
#endif
        }
}
