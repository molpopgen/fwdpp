/*! \include diploid_ind.cc
  Simulate a single, finite Wright-Fisher diploid population with mutation,
  recombination, and no selection.

  This program illustrates many features of fwdpp:
  1.  Custom mutation classes
  2.  Implementing a mutation model (infinitely-many sites)
  3.  Iterating a population through its life cycle
  4.  Outputting a sample in "ms" format
*/
#include <iostream>
#include <unordered_map>
#include <type_traits>
#include <vector>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include <fwdpp/diploid.hh>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
// typedef mutation_with_age mtype;
using mtype = fwdpp::popgenmut;
#define SINGLEPOP_SIM
#include <common_ind.hpp>
#include <numeric>
void
compact(singlepop_t &pop)
{
    std::vector<fwdpp::uint_t> indexes(pop.mutations.size());
    std::iota(std::begin(indexes), std::end(indexes), 0);

    auto new_indexes_end = std::stable_partition(
        std::begin(indexes), std::end(indexes),
        [&pop](const fwdpp::uint_t i) { return pop.mcounts[i]; });

    std::sort(std::begin(indexes), new_indexes_end,
              [&pop](const fwdpp::uint_t i, const fwdpp::uint_t j) {
                  return pop.mutations[i].pos < pop.mutations[j].pos;
              });
    std::vector<fwdpp::uint_t> reindex(indexes.size());
    std::size_t new_indexes_size
        = std::distance(std::begin(indexes), new_indexes_end);
    for (std::size_t i = 0; i < new_indexes_size; ++i)
        {
            reindex[indexes[i]] = i;
        }
    for (auto &g : pop.gametes)
        {
            if (g.n)
                {
                    for (auto &m : g.mutations)
                        {
                            m = reindex[m];
                        }
                    for (auto &m : g.smutations)
                        {
                            m = reindex[m];
                        }
                }
        }
    decltype(pop.mutations) reordered_muts;
    decltype(pop.mcounts) reordered_mcounts;
    reordered_muts.reserve(pop.mutations.size());
    reordered_mcounts.reserve(pop.mutations.size());
    for (auto i : indexes)
        {
            reordered_muts.emplace_back(std::move(pop.mutations[i]));
            reordered_mcounts.push_back(pop.mcounts[i]);
            if (reordered_mcounts.back() > 0)
                {
                    auto x = pop.mut_lookup.equal_range(
                        reordered_muts.back().pos);
                    while (x.first != x.second)
                        {
                            if (x.first->second == i)
                                {
                                    x.first->second
                                        = reordered_muts.size() - 1;
                                    break;
                                }
                            ++x.first;
                        }
                }
        }
    pop.mutations.swap(reordered_muts);
    pop.mcounts.swap(reordered_mcounts);
}

int
main(int argc, char **argv)
{
    if (argc != 8)
        {
            std::cerr << "Too few arguments\n"
                      << "Usage: diploid_ind N theta rho ngens samplesize "
                         "nreps seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N = unsigned(atoi(argv[argument++])); // Number of diploids
    const double theta = atof(argv[argument++]); // 4*n*mutation rate.  Note:
                                                 // mutation rate is per
                                                 // REGION, not SITE!!
    const double rho = atof(argv[argument++]);   // 4*n*recombination rate.
                                                 // Note: recombination rate is
                                                 // per REGION, not SITE!!
    const unsigned ngens = unsigned(
        atoi(argv[argument++])); // Number of generations to simulate
    const unsigned samplesize1 = unsigned(
        atoi(argv[argument++])); // Sample size to draw from the population
    int nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned seed
        = unsigned(atoi(argv[argument++])); // Random number seed

    const double mu = theta / double(4 * N); // per-gamete mutation rate

    /*
      littler r is the recombination rate per region per generation.
    */
    const double littler = rho / double(4 * N);

    // Write the command line to stderr
    std::copy(argv, argv + argc,
              std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    // Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);

    unsigned twoN = 2 * N;

    // recombination map is uniform[0,1)
    const auto rec
        = fwdpp::recbinder(fwdpp::poisson_xover(littler, 0., 1.), r.get());

    while (nreps--)
        {
            singlepop_t pop(N);
            pop.mutations.reserve(
                size_t(std::ceil(std::log(2 * N) * theta + 0.667 * theta)));
            unsigned generation = 0;
            double wbar;

            const auto mmodel =
                [&pop, &r, &generation](std::queue<std::size_t> &recbin,
                                        singlepop_t::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, r.get(), pop.mut_lookup, generation,
                        0.0, [&r]() { return gsl_rng_uniform(r.get()); },
                        []() { return 0.0; }, []() { return 0.0; });
                };
            for (generation = 0; generation < ngens; ++generation)
                {
                    // Iterate the population through 1 generation
                    wbar = fwdpp::sample_diploid(
                        r.get(),
                        pop.gametes,   // non-const reference to gametes
                        pop.diploids,  // non-const reference to diploids
                        pop.mutations, // non-const reference to mutations
                        pop.mcounts,
                        N,  // current pop size, remains constant
                        mu, // mutation rate per gamete
                        /*
                          The mutation model will be applied
                          by
                          sample_diploid in order to add mutations to gametes
                          each generation.
                        */
                        mmodel,
                        // The function to generation recombination positions:
                        rec,
                        /*
                          Fitness is multiplicative over sites.

                          The fitness model takes two gametes and a
                          vector of mutations as arguments.
                          The gametes are passed to this function by
                          fwdpp::sample_diploid, and the _1 and _2 are
                          placeholders for
                          those gametes (see documentation for boost/bind.hpp
                          for details).
                          The mutation container is passed in as _3.
                          The 2. means that fitnesses are 1, 1+sh, and 1+2s for
                          genotypes
                          AA, Aa, and aa, respectively, where a is a mutation
                          with
                          selection coefficient s and dominance h, and the
                          fitness of
                          the diploid is the product of fitness over sites

                          There is no selection in this simulation, but this
                          function is called anyways to illustrate it as
                          multiplicative
                          models are very common in population genetics
                        */
                        fwdpp::multiplicative_diploid(), pop.neutral,
                        pop.selected);
                    fwdpp::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, twoN);
                    if (generation && generation % 100 == 0.0)
                        {
                            compact(pop);
                        }
                    assert(fwdpp::check_sum(pop.gametes, twoN));
                    //for(std::size_t i=0;i<pop.mcounts.size();++i)
                    //{
                    //	if(pop.mcounts[i])
                    //	{
                    //		if(pop.mut_lookup.find(pop.mutations[i].pos)==pop.mut_lookup.end())
                    //		{
                    //			throw 1;
                    //		}
                    //	}
                    //}
                }

            // Take a sample of size samplesize1 from the population
            std::vector<std::pair<double, std::string>> mslike
                = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                                   pop.diploids, samplesize1, true);

// Write the sample date a to libsequence's Sequence::SimData and
// print to screen
#ifdef HAVE_LIBSEQUENCE
            Sequence::SimData sdata;
            if (!mslike.empty())
                {
                    sdata.assign(mslike.begin(), mslike.end());
                    std::cout << sdata << '\n';
                }
            else
                {
                    std::cout << "//\nsegsites: 0\n";
                }
#endif
        }
    return 0;
}
