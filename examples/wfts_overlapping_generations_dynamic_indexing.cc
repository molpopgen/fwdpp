#include <iostream>
#include <cstdlib>
#include "wfevolvets.hpp"

const double genome_length = 1.0;

int
main(int argc, char** argv)
{
    if (argc < 8)
        {
            std::cerr << "usage: " << argv[0]
                      << "popsize ngenerations psurvival recrate "
                         "check_interval simplify seed\n";
            std::exit(1);
        }
    int argnext = 1;
    unsigned N = std::atoi(argv[argnext++]);
    unsigned ngenerations = std::atoi(argv[argnext++]);
    double psurvival = std::atof(argv[argnext++]);
    double recrate = std::atof(argv[argnext++]);
    unsigned check_interval = std::atoi(argv[argnext++]);
    unsigned simplify = std::atoi(argv[argnext++]);
    unsigned seed = std::atoi(argv[argnext++]);

    fwdpp::poisson_interval recombination(0, genome_length, recrate);
    GSLrng rng(seed);
    // Initialize a table_collection with 2N nodes at time zero,
    // in pop zero, and with our genome length
    fwdpp::ts::table_collection tables(2 * N, 0, 0, genome_length);

    // Set up individual metadata
    std::vector<diploid_metadata> metadata;
    for (unsigned i = 0; i < N; ++i)
        {
            // label, time, fitness, node 1, node 2
            metadata.emplace_back(i, 0, 1.0, 2 * i, 2 * i + 1);
        }

    wfevolvets_dynamic_indexing(rng, ngenerations, 0, false, false, check_interval, simplify,
                                psurvival, 0.0, recombination, metadata,
                                tables);
}
