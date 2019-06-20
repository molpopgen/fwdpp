/*  \include juvenile_migration.cc
 *
 *  Juvenile migration.
 *  Selection.
 *  Constant pop size.
 *  Unequal migration rates.
 *
 * The goal is to show how to use a population
 * object to handle details of population structure.
 *
 * The outline of the scheme is as follows:
 *
 * 1. Offspring are generated in sorted order by deme
 * label.  In other words, deme 1 before deme 2.
 *
 * 2. Each generation, we migrate first.  Thus,
 * we may start a generation with parents like this:
 *
 * Parent Deme
 * 0    0
 * 1    0
 * 2    0
 * 3    1
 * 4    1
 * 5    1
 *
 * After migration, we may end up with:
 *
 * Parent Deme
 * 0    1*
 * 1    0
 * 2    0
 * 3    1
 * 4    0*
 * 5    1
 *
 * An asterisk represent migrant individuals.
 *
 * We will generate efficient lookup tables to sample individuals
 * proportional to their within-deme fitness, post-migration.
 *
 * These lookup tables will map to the indexes of the parents of each
 * deme:
 *
 * Lookup 1:
 *
 * Index Parent
 * 0    1
 * 1    2
 * 2    4
 *
 * Lookup 2:
 * Index Parent
 * 0    0
 * 1    3
 * 1    5
 *
 * The code is not safe for real-world use.  It doesn't handle the corner
 * case of all of deme 1 migrating into deme 2, leaving deme 1 "extinct"
 * in the next generation.  Such edge cases are clearly important, but
 * the goal here it to show the book-keeping of parental fitnesses
 * and the mapping back to deme labels.
 */
#include <fwdpp/diploid.hh>
#include <fwdpp/recbinder.hpp>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include <numeric>
#include <functional>
#include <cassert>
#include <unordered_set>
#include <fwdpp/debug.hpp>
#include <fwdpp/popgenmut.hpp>
#define DIPLOID_POPULATION_SIM
// the type of mutation
using mtype = fwdpp::popgenmut;
#include <common_ind.hpp>
#include <gsl/gsl_randist.h>

struct parent_lookup_tables
// Our object for choosing parents each generation
{
    // These return indexes of parents from demes 1 and 2,
    // resp, chosen in O(1) time proportional to
    // relative fitness within each deme
    fwdpp::gsl_ran_discrete_t_ptr lookup1, lookup2;
    // These vectors map indexes returned from sampling
    // lookup1 and lookup2 to diploids in the population
    // object.
    std::vector<std::size_t> parents1, parents2;
};

template <typename fitness_fxn>
parent_lookup_tables
migrate_and_calc_fitness(const gsl_rng *r, diploid_population &pop,
                         const fitness_fxn &wfxn, const fwdpp::uint_t N1,
                         const fwdpp::uint_t N2, const double m12,
                         const double m21)
// This function will be called at the start of each generation.
// The main goal is to return the lookup tables described above.
// But, "while we're at it", it does some other stuff that
// needs to be done at the start of each generation.
// Neither the most rigorous nor the most efficient:
// 1. Ignores probability of back-migration.
// 2. Allocates 4 vectors each generation.
{
    parent_lookup_tables rv;

    // Temp containers for fitnesses in each deme,
    // post-migration
    std::vector<double> w1, w2;

    // Pick no. migrants 1 -> 2 and 2 -> 1.
    unsigned nmig12 = gsl_ran_poisson(r, static_cast<double>(N1) * m12);
    unsigned nmig21 = gsl_ran_poisson(r, static_cast<double>(N2) * m21);

    // Fill a vector of N1 zeros and N2 ones:
    std::vector<fwdpp::uint_t> deme_labels(N1, 0);
    deme_labels.resize(N1 + N2, 1);
    assert(deme_labels.size() == pop.diploids.size());

    // Set up source and destination containers
    // for sampling w/o replacement
    std::vector<std::size_t> individuals(N1 + N2);
    std::iota(std::begin(individuals), std::end(individuals), 0);
    std::vector<std::size_t> migrants(std::max(nmig12, nmig21));

    // Who is migrating 1 -> 2?
    gsl_ran_choose(r, migrants.data(), nmig12, individuals.data(), N1,
                   sizeof(std::size_t));
    for (std::size_t i = 0; i < nmig12; ++i)
        {
            deme_labels[migrants[i]] = !deme_labels[migrants[i]];
        }

    // Exact same logic for migrants 2 -> 1
    gsl_ran_choose(r, migrants.data(), nmig21, individuals.data() + N1, N2,
                   sizeof(std::size_t));
    for (std::size_t i = 0; i < nmig21; ++i)
        {
            deme_labels[migrants[i]] = !deme_labels[migrants[i]];
        }

    // Go over all parents, set haploid_genomes counts to zero,
    // and put individual IDs and fitnesses into
    // the right vectors:
    for (std::size_t i = 0; i < deme_labels.size(); ++i)
        {
            // fwdpp requires that we zero out haploid_genome
            // counts each generation.  Since we're looping
            // over diploids here, now is a good time to
            // handle this task, which saves us from having to
            // do another O(N1+N2) loop:
            pop.haploid_genomes[pop.diploids[i].first].n
                = pop.haploid_genomes[pop.diploids[i].second].n = 0;
            if (deme_labels[i] == 0)
                {
                    rv.parents1.push_back(i);
                    w1.push_back(
                        wfxn(pop.diploids[i], pop.haploid_genomes, pop.mutations));
                }
            else
                {
                    rv.parents2.push_back(i);
                    w2.push_back(
                        wfxn(pop.diploids[i], pop.haploid_genomes, pop.mutations));
                }
        }

    // Set up our lookup tables:
    rv.lookup1.reset(gsl_ran_discrete_preproc(rv.parents1.size(), w1.data()));
    rv.lookup2.reset(gsl_ran_discrete_preproc(rv.parents2.size(), w2.data()));
    return rv;
};

template <typename fitness_fxn, typename rec_fxn, typename mut_fxn>
void
evolve_two_demes(const gsl_rng *r, diploid_population &pop,
                 const fwdpp::uint_t N1, const fwdpp::uint_t N2,
                 const double m12, const double m21, const double mu,
                 const fitness_fxn &wfxn, const rec_fxn &recfxn,
                 const mut_fxn &mutfxn)
{
    // Handle mutation/haploid_genome "recycling":
    auto mut_recycling_bin = fwdpp::make_mut_queue(pop.mcounts);
    auto gam_recycling_bin = fwdpp::make_haploid_genome_queue(pop.haploid_genomes);

    // Migration and build lookup tables:
    auto lookups = migrate_and_calc_fitness(r, pop, wfxn, N1, N2, m12, m21);

#ifndef NDEBUG
    for (const auto &g : pop.haploid_genomes)
        assert(!g.n);
#endif

    // Copy parents
    const auto parents(pop.diploids);

    // Fill in the next generation!
    // We generate the offspring for deme 1 first,
    // and then for deme 2
    for (fwdpp::uint_t i = 0; i < N1 + N2; ++i)
        {
            std::size_t p1 = std::numeric_limits<std::size_t>::max();
            std::size_t p2 = std::numeric_limits<std::size_t>::max();
            if (i < N1) // pick parents from pop 1
                {
                    p1 = lookups.parents1[gsl_ran_discrete(
                        r, lookups.lookup1.get())];
                    p2 = lookups.parents1[gsl_ran_discrete(
                        r, lookups.lookup1.get())];
                }
            else // pick parents from pop 2
                {
                    p1 = lookups.parents2[gsl_ran_discrete(
                        r, lookups.lookup2.get())];
                    p2 = lookups.parents2[gsl_ran_discrete(
                        r, lookups.lookup2.get())];
                }
            assert(p1 < parents.size());
            assert(p2 < parents.size());

            /*
              These are the haploid_genomes from each parent.
            */
            auto p1g1 = parents[p1].first;
            auto p1g2 = parents[p1].second;
            auto p2g1 = parents[p2].first;
            auto p2g2 = parents[p2].second;

            // "Mendel"
            if (gsl_rng_uniform(r) < 0.5)
                std::swap(p1g1, p1g2);
            if (gsl_rng_uniform(r) < 0.5)
                std::swap(p2g1, p2g2);

            mutate_recombine_update(r, pop.haploid_genomes, pop.mutations,
                                    std::make_tuple(p1g1, p1g2, p2g1, p2g2),
                                    recfxn, mutfxn, mu, gam_recycling_bin,
                                    mut_recycling_bin, pop.diploids[i],
                                    pop.neutral, pop.selected);
        }
    fwdpp::debug::validate_sum_haploid_genome_counts(pop.haploid_genomes,
                                             2 * pop.diploids.size());
#ifndef NDEBUG
    for (const auto &dip : pop.diploids)
        {
            assert(pop.haploid_genomes[dip.first].n > 0);
            assert(pop.haploid_genomes[dip.first].n <= 2 * (N1 + N2));
            assert(pop.haploid_genomes[dip.second].n > 0);
            assert(pop.haploid_genomes[dip.second].n <= 2 * (N1 + N2));
        }
#endif

    // Update mutation counts
    fwdpp::fwdpp_internal::process_haploid_genomes(pop.haploid_genomes, pop.mutations,
                                           pop.mcounts);

    assert(pop.mcounts.size() == pop.mutations.size());
#ifndef NDEBUG
    for (const auto &mc : pop.mcounts)
        {
            assert(mc <= 2 * (N1 + N2));
        }
#endif
    fwdpp::debug::validate_pop_data(pop);

    // Prune fixations from haploid_genomes
    fwdpp::fwdpp_internal::haploid_genome_cleaner(pop.haploid_genomes, pop.mutations,
                                          pop.mcounts, 2 * (N1 + N2),
                                          std::true_type());
}

int
main(int argc, char **argv)
{
    if (argc != 12)
        {
            std::cerr
                << "Too few arguments.\n"
                << "Usage: juvenile_migration N1 N1 m12 m21 theta_neutral "
                   "theta_deleterious rho s h ngens "
                   "seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N1 = atoi(argv[argument++]);
    const unsigned N2 = atoi(argv[argument++]);
    const double m12 = atof(argv[argument++]);
    const double m21 = atof(argv[argument++]);
    const double theta_neutral = atof(argv[argument++]);
    const double theta_del = atof(argv[argument++]);
    const double rho = atof(argv[argument++]);
    const double s = atof(argv[argument++]);
    const double h = atof(argv[argument++]);
    const unsigned ngens = atoi(argv[argument++]);
    const unsigned seed = atoi(argv[argument++]);

    const unsigned N = N1 + N2; // Total metapop size, for convenience
    const double mu_neutral = theta_neutral / double(4 * N);
    const double mu_del = theta_del / double(4 * N);
    const double littler = rho / double(4 * N);

    std::copy(argv, argv + argc,
              std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    GSLrng r(seed);

    // recombination map is uniform[0,1)
    const auto rec
        = fwdpp::recbinder(fwdpp::poisson_xover(littler, 0., 1.), r.get());

    const double pselected = mu_del / (mu_del + mu_neutral);

    auto wfxn = fwdpp::multiplicative_diploid(fwdpp::fitness(1.));
    diploid_population pop(N);
    pop.mutations.reserve(
        size_t(std::ceil(std::log(2 * N) * (theta_neutral + theta_del)
                         + 0.667 * (theta_neutral + theta_del))));
    unsigned generation = 0;
    const auto mmodel = [&pop, &r, &generation, s, h,
                         pselected](fwdpp::flagged_mutation_queue &recbin,
                                    diploid_population::mcont_t &mutations) {
        return fwdpp::infsites_popgenmut(
            recbin, mutations, r.get(), pop.mut_lookup, generation, pselected,
            [&r]() { return gsl_rng_uniform(r.get()); }, [s]() { return s; },
            [h]() { return h; });
    };

    double wbar = 1;
    for (generation = 0; generation < ngens; ++generation)
        {
            fwdpp::debug::validate_sum_haploid_genome_counts(pop.haploid_genomes,
                                                     2 * (N1 + N2));

            // Call our fancy new evolve function
            evolve_two_demes(r.get(), pop, N1, N2, m12, m21,
                             mu_neutral + mu_del, wfxn, rec, mmodel);
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
            fwdpp::debug::validate_sum_haploid_genome_counts(pop.haploid_genomes, 2 * N);
        }
}
