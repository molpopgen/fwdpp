/*
  \include diploid_gzbinaryIO_ind.cc

  Same as diploid_binaryIO_ind.cc, but input/output are gzip compressed

  The population is then read back in and compared to what was written out.

  Main point here is to show how the write/read function objects for
  fwdpp/IO.hpp should be written.

  Also illustrates POSIX file locking via <fcntl.h>, which is super-useful on
  clusters.

  Example use that runs quickly: ./diploid_gzbinaryIO 1000 10 10 10000 1 index
  haps.bin.gz $RANDOM
*/

#include <fwdpp/diploid.hh>
#include <numeric>
#include <functional>
#include <cassert>
#include <fstream>
#include <fcntl.h>
#include <fwdpp/sugar/infsites.hpp>
using mtype = KTfwd::popgenmut;
#define SINGLEPOP_SIM
#include <common_ind.hpp>

int
main(int argc, char **argv)
{
    if (argc != 9)
        {
            std::cerr << "Too few arguments\n"
                      << "Usage: diploid_gzbinaryIO_ind N theta rho ngens "
                         "replicate_no indexfile hapfile seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N = atoi(argv[argument++]);
    const double theta = atof(argv[argument++]);
    const double rho = atof(argv[argument++]);
    const unsigned ngens = atoi(argv[argument++]);
    const unsigned replicate_no = atoi(argv[argument++]);
    const char *indexfile = argv[argument++];
    const char *hapfile = argv[argument++];
    const unsigned seed = atoi(argv[argument++]);

    const double mu = theta / double(4 * N);
    const double littler = rho / double(4 * N);

    std::copy(argv, argv + argc,
              std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    GSLrng r(seed);
    singlepop_t pop(N);
    pop.mcounts.reserve(
        size_t(std::ceil(std::log(2 * N) * theta + 0.667 * theta)));
    unsigned generation;
    double wbar;

    // recombination map is uniform[0,1)
    std::function<double(void)> recmap = std::bind(gsl_rng_uniform, r.get());

    for (generation = 0; generation < ngens; ++generation)
        {
            wbar = KTfwd::sample_diploid(
                r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                N, mu, std::bind(KTfwd::infsites(), std::placeholders::_1,
                                 std::placeholders::_2, r.get(),
                                 std::ref(pop.mut_lookup), generation, mu, 0.,
                                 [&r]() { return gsl_rng_uniform(r.get()); },
                                 []() { return 0.; }, []() { return 0.; }),
                // The function to generation recombination positions:
                std::bind(KTfwd::poisson_xover(), r.get(), littler, 0., 1.,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                std::bind(KTfwd::multiplicative_diploid(),
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, 2.),
                pop.neutral, pop.selected);
            KTfwd::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
        }

    /*
      Note: gzFiles are not as easily compatible with file-locking and creating
      an index, etc.,
      as done in diploid_binaryIO.cc.  See
      https://github.com/molpopgen/BigDataFormats
      for why not.
     */

    // establish POSIX file locks for output
    struct flock index_flock;
    index_flock.l_type = F_WRLCK; /*Write lock*/
    index_flock.l_whence = SEEK_SET;
    index_flock.l_start = 0;
    index_flock.l_len = 0; /*Lock whole file*/

    FILE *index_fh = fopen(indexfile, "a");
    int index_fd = fileno(index_fh);
    if (index_fd == -1)
        {
            std::cerr << "ERROR: could not open " << indexfile << '\n';
            exit(0);
        }
    if (fcntl(index_fd, F_SETLKW, &index_flock) == -1)
        {
            std::cerr << "ERROR: could not obtain lock on " << indexfile
                      << '\n';
            exit(0);
        }

    // OK, we no have an exclusive lock on the index file.  In principle, that
    // is sufficient for us to move on.

    // Now, write output to a gzipped file
    gzFile gzout = gzopen(hapfile, "ab"); // open in append mode.  The b =
                                          // binary mode.  Not required on all
                                          // systems, but never hurts.
    long long written
        = KTfwd::gzserialize()(gzout, pop, KTfwd::mutation_writer());
    // write info to index file
    fprintf(index_fh, "%u %ld\n", replicate_no, gztell(gzout));
    gzclose(gzout);

    // We can now release close the index file, release the lock, etc.
    index_flock.l_type = F_UNLCK;
    if (fcntl(index_fd, F_UNLCK, &index_flock) == -1)
        {
            std::cerr << "ERROR: could not release lock on " << indexfile
                      << '\n';
            exit(0);
        }
    fflush(index_fh);
    fclose(index_fh);

    /*
     Read the data back in.
     Unlike our other examples, we must keep track of the cumulative offsets in
     the index file
     up until we find our record_number
   */
    std::ifstream index_in(indexfile);
    unsigned long long rec_offset = 0, offset_i;
    unsigned rid;
    while (!index_in.eof())
        {
            index_in >> rid >> offset_i >> std::ws;
            if (rid != replicate_no)
                {
                    rec_offset += offset_i;
                }
            else
                {
                    break;
                }
        }
    index_in.close();
    if (rid != replicate_no)
        {
            std::cerr << "Error: replicate id not found in index file, "
                      << replicate_no << ' ' << rid << '\n';
            exit(0);
        }

    singlepop_t pop2(pop.N);

    gzFile gzin
        = gzopen(hapfile, "rb"); // open it for reading.  Again, binary mode.
    gzseek(gzin, rec_offset, SEEK_CUR); // seek to position

    KTfwd::gzdeserialize()(
        pop2, gzin,
        std::bind(KTfwd::mutation_reader<mtype>(), std::placeholders::_1));

    for (std::size_t i = 0; i < pop.diploids.size(); ++i)
        {
            std::cout << "Diploid no. " << i << ": " << pop.diploids[i].first
                      << ',' << pop.diploids[i].second << ' '
                      << pop2.diploids[i].first << ','
                      << pop2.diploids[i].second << " -> "
                      << (pop.gametes[pop.diploids[i].first].mutations
                          == pop2.gametes[pop2.diploids[i].first].mutations)
                      << ' '
                      << (pop.gametes[pop.diploids[i].first].smutations
                          == pop2.gametes[pop2.diploids[i].first].smutations)
                      << ' '
                      << (pop.gametes[pop.diploids[i].second].mutations
                          == pop2.gametes[pop2.diploids[i].second].mutations)
                      << ' '
                      << (pop.gametes[pop.diploids[i].second].smutations
                          == pop2.gametes[pop2.diploids[i].second].smutations)
                      << '\n';
        }
    std::cout << "Are mutation containers the same: "
              << (pop.mutations == pop2.mutations) << '\n';

    // Now, compare what we wrote to what we read
    // std::cout << pop.gametes.size() << ' ' << pop2.gametes.size() << ' ' <<
    // pop.mutations.size() << ' ' << pop2.mutations.size()
    // 	    << ' ' << pop.diploids.size() << ' ' << pop2.diploids.size() <<
    // '\n';

    // for( unsigned i = 0 ; i < pop.diploids.size() ; ++i )
    //   {
    //     std::cout << "Diploid " << i << ":\nWritten:\n";
    //     for( unsigned j = 0 ; j < pop.diploids[i].first->mutations.size() ;
    //     ++j )
    // 	{
    // 	  std::cout << '(' << pop.diploids[i].first->mutations[j]->pos << ','
    // 		    << pop.diploids[i].first->mutations[j]->n << ')';
    // 	}
    //     std::cout << '\n';
    //     for( unsigned j = 0 ; j < pop.diploids[i].second->mutations.size() ;
    //     ++j )
    // 	{
    // 	  std::cout << '(' << pop.diploids[i].second->mutations[j]->pos << ','
    // 		    << pop.diploids[i].second->mutations[j]->n << ')';
    // 	}
    //     std::cout << '\n';

    //     std::cout << "Read:\n";
    //     for( unsigned j = 0 ; j < pop2.diploids[i].first->mutations.size() ;
    //     ++j )
    // 	{
    // 	  std::cout << '(' << pop2.diploids[i].first->mutations[j]->pos << ','
    // 		    << pop2.diploids[i].first->mutations[j]->n << ')';
    // 	}
    //     std::cout << '\n';
    //     for( unsigned j = 0 ; j < pop2.diploids[i].second->mutations.size()
    //     ; ++j )
    // 	{
    // 	  std::cout << '(' << pop2.diploids[i].second->mutations[j]->pos << ','
    // 		    << pop2.diploids[i].second->mutations[j]->n << ')';
    // 	}
    //     std::cout << '\n';
    //   }
}
