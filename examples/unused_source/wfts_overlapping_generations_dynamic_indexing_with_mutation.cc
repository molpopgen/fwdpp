#include <iostream>
#include <cstdlib>
#include "wfevolvets.hpp"

const double genome_length = 1.0;

int
main(int argc, char** argv)
{
    if (argc < 11)
        {
            std::cerr << "usage: " << argv[0]
                      << " popsize ngenerations psurvival recrate mutrate "
                         "count_mutations_interval remove_fixations "
                         "track_samples simplify seed\n";
            std::exit(1);
        }
    int argnext = 1;
    unsigned N = std::atoi(argv[argnext++]);
    unsigned ngenerations = std::atoi(argv[argnext++]);
    double psurvival = std::atof(argv[argnext++]);
    double recrate = std::atof(argv[argnext++]);
    double mutrate = std::atof(argv[argnext++]);
    unsigned count_mutations_interval = std::atoi(argv[argnext++]);
    bool remove_fixations = std::atoi(argv[argnext++]);
    bool track_samples = std::atoi(argv[argnext++]);
    unsigned simplify = std::atoi(argv[argnext++]);
    unsigned seed = std::atoi(argv[argnext++]);

    fwdpp::poisson_interval recombination(0, genome_length, recrate);
    GSLrng rng(seed);
    // Initialize a table_collection with 2N nodes at time zero,
    // in pop zero, and with our genome length
    fwdpp::ts::std_table_collection tables(2 * N, 0, 0, genome_length);

    // Set up individual metadata
    std::vector<diploid_metadata> metadata;
    for (unsigned i = 0; i < N; ++i)
        {
            // label, time, fitness, node 1, node 2
            metadata.emplace_back(i, 0, 1.0, 2 * i, 2 * i + 1);
        }

    auto mcounts = wfevolvets_dynamic_indexing(
        rng, ngenerations, count_mutations_interval, remove_fixations,
        track_samples, 0, simplify, psurvival, mutrate, recombination,
        metadata, tables);
    if (!mcounts.empty())
        {
            auto itr = std::max_element(begin(mcounts), end(mcounts));
            decltype(mcounts) afs((*itr) + 1, 0);
            for (auto m : mcounts)
                {
                    if (m >= afs.size())
                        {
                            throw std::runtime_error("bad afs index");
                        }
                    afs[m]++;
                }
            for (std::size_t i = 0; i < afs.size(); ++i)
                {
                    if (afs[i] > 0)
                        {
                            std::cout << i << ' ' << afs[i] << '\n';
                        }
                }
        }
}

