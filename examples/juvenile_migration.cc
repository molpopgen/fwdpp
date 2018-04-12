/*
  \include juvenile_migration.cc
 *  Juvenile migration.
 *  Selection.
 *  Constant pop size.
*/
#include <fwdpp/diploid.hh>
#include <fwdpp/recbinder.hpp>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include <numeric>
#include <functional>
#include <cassert>
#include <iomanip>
#include <fwdpp/sugar/popgenmut.hpp>
#define SINGLEPOP_SIM
// the type of mutation
using mtype = fwdpp::popgenmut;
#include <common_ind.hpp>
#include <gsl/gsl_randist.h>
using namespace fwdpp;

struct parent_lookup_tables
{
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup1, lookup2;
    std::vector<std::size_t> parents1, parents2;
};

template <typename fitness_fxn>
parent_lookup_tables
migrate_and_calc_fitness(const gsl_rng *r, singlepop_t &pop,
                         const fitness_fxn &wfxn, const uint_t N1,
                         const uint_t N2, const double m12, const double m21)
// Neither the most rigorous nor the most efficient:
// 1. Ignores probability of back-migration.
// 2. Allocates 4 vectors each generation.
{
    parent_lookup_tables rv;
    std::vector<double> w1, w2;

    unsigned nmig12 = gsl_ran_poisson(r, static_cast<double>(N1) * m12);
    unsigned nmig21 = gsl_ran_poisson(r, static_cast<double>(N2) * m21);

    std::vector<uint_t> deme_labels(N1, 0);
    deme_labels.resize(N1 + N2, 1);
    assert(deme_labels.size() == pop.diploids.size());
    std::unordered_set<uint_t> migrants;
    for (unsigned i = 0; i < nmig12; ++i)
        {
            auto mig = static_cast<uint_t>(gsl_ran_flat(r, 0, N1));
            while (migrants.find(mig) != migrants.end())
                {
                    mig = static_cast<uint_t>(gsl_ran_flat(r, 0, N1));
                }
            deme_labels[i] = !deme_labels[i];
            migrants.insert(mig);
        }
    for (unsigned i = 0; i < nmig21; ++i)
        {
            auto mig = static_cast<uint_t>(gsl_ran_flat(r, N1, N1 + N2));
            while (migrants.find(mig) != migrants.end())
                {
                    mig = static_cast<uint_t>(gsl_ran_flat(r, N1, N1 + N2));
                }
            deme_labels[i] = !deme_labels[i];
            migrants.insert(mig);
        }
    for (std::size_t i = 0; i < deme_labels.size(); ++i)
        {
            pop.gametes[pop.diploids[i].first].n
                = pop.gametes[pop.diploids[i].second].n = 0;
            if (deme_labels[i] == 0)
                {
                    rv.parents1.push_back(i);
                    w1.push_back(
                        wfxn(pop.diploids[i], pop.gametes, pop.mutations));
                }
            else
                {
                    rv.parents2.push_back(i);
                    w2.push_back(
                        wfxn(pop.diploids[i], pop.gametes, pop.mutations));
                }
        }
    rv.lookup1.reset(gsl_ran_discrete_preproc(rv.parents1.size(), w1.data()));
    rv.lookup2.reset(gsl_ran_discrete_preproc(rv.parents2.size(), w2.data()));
    return rv;
};

template <typename fitness_fxn, typename rec_fxn, typename mut_fxn>
void
evolve_two_demes(const gsl_rng *r, singlepop_t &pop, const uint_t N1,
                 const uint_t N2, const double m12, const double m21,
                 const double mu, const fitness_fxn &wfxn,
                 const rec_fxn &recfxn, const mut_fxn &mutfxn)
{
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(pop.mcounts);
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(pop.gametes);

    auto lookups = migrate_and_calc_fitness(r, pop, wfxn, N1, N2, m12, m21);
#ifndef NDEBUG
    for (const auto &g : pop.gametes)
        assert(!g.n);
#endif
    const auto parents(pop.diploids);

    // Fill in the next generation!
    for (uint_t i = 0; i < N1 + N2; ++i)
        {
            std::size_t p1 = std::numeric_limits<std::size_t>::max();
            std::size_t p2 = std::numeric_limits<std::size_t>::max();
            if (i < N1)
                {
                    p1 = lookups.parents1[gsl_ran_discrete(
                        r, lookups.lookup1.get())];
                    p2 = lookups.parents1[gsl_ran_discrete(
                        r, lookups.lookup1.get())];
                }
            else
                {
                    p1 = lookups.parents2[gsl_ran_discrete(
                        r, lookups.lookup2.get())];
                    p2 = lookups.parents2[gsl_ran_discrete(
                        r, lookups.lookup2.get())];
                }
            assert(p1 < parents.size());
            assert(p2 < parents.size());
            /*
              These are the gametes from each parent.
              This is a trivial assignment if keys.
            */
            auto p1g1 = parents[p1].first;
            auto p1g2 = parents[p1].second;
            auto p2g1 = parents[p2].first;
            auto p2g2 = parents[p2].second;

            if (gsl_rng_uniform(r) < 0.5)
                std::swap(p1g1, p1g2);
            if (gsl_rng_uniform(r) < 0.5)
                std::swap(p2g1, p2g2);

            mutate_recombine_update(r, pop.gametes, pop.mutations,
                                    std::make_tuple(p1g1, p1g2, p2g1, p2g2),
                                    recfxn, mutfxn, mu, gam_recycling_bin,
                                    mut_recycling_bin, pop.diploids[i],
                                    pop.neutral, pop.selected);
        }
    assert(check_sum(pop.gametes, 2 * (N1 + N2)));
#ifndef NDEBUG
    for (const auto &dip : pop.diploids)
        {
            assert(pop.gametes[dip.first].n > 0);
            assert(pop.gametes[dip.first].n <= 2 * (N1 + N2));
            assert(pop.gametes[dip.second].n > 0);
            assert(pop.gametes[dip.second].n <= 2 * (N1 + N2));
        }
#endif
    fwdpp_internal::process_gametes(pop.gametes, pop.mutations, pop.mcounts);

    assert(pop.mcounts.size() == pop.mutations.size());
#ifndef NDEBUG
    for (const auto &mc : pop.mcounts)
        {
            assert(mc <= 2 * (N1 + N2));
        }
#endif
    assert(
        popdata_sane(pop.diploids, pop.gametes, pop.mutations, pop.mcounts));

    fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts,
                                   2 * (N1 + N2), std::true_type());
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

    auto wfxn = fwdpp::multiplicative_diploid(1.);
    singlepop_t pop(N);
    pop.mutations.reserve(
        size_t(std::ceil(std::log(2 * N) * (theta_neutral + theta_del)
                         + 0.667 * (theta_neutral + theta_del))));
    unsigned generation = 0;
    const auto mmodel = [&pop, &r, &generation, s, h,
                         pselected](std::queue<std::size_t> &recbin,
                                    singlepop_t::mcont_t &mutations) {
        return fwdpp::infsites_popgenmut(
            recbin, mutations, r.get(), pop.mut_lookup, generation, pselected,
            [&r]() { return gsl_rng_uniform(r.get()); }, [s]() { return s; },
            [h]() { return h; });
    };

    double wbar = 1;
    for (generation = 0; generation < ngens; ++generation)
        {
            assert(fwdpp::check_sum(pop.gametes, 2 * (N1 + N2)));
            evolve_two_demes(r.get(), pop, N1, N2, m12, m21,
                             mu_neutral + mu_del, wfxn, rec, mmodel);
            wbar = fwdpp::sample_diploid(
                r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                N, mu_neutral + mu_del, mmodel,
                // The function to generation recombination positions:
                rec, wfxn, pop.neutral, pop.selected);
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation, 2 * N);
            assert(fwdpp::check_sum(pop.gametes, 2 * N));
        }
}
