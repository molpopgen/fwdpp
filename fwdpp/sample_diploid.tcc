//  -*- C++ -*-

/*!
  \file sample_diploid.tcc

  \brief Definitions of functions for evolving populations of diploids.
*/

#ifndef FWDPP_SAMPLE_DIPLOID_TCC
#define FWDPP_SAMPLE_DIPLOID_TCC

#include <cassert>
#include <fwdpp/debug.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace fwdpp
{
    // single deme, constant N
    template <typename gamete_type, typename gamete_cont_type_allocator,
              typename mutation_type, typename mutation_cont_type_allocator,
              typename diploid_geno_t, typename diploid_vector_type_allocator,
              typename diploid_fitness_function, typename mutation_model,
              typename recombination_policy,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              typename mutation_removal_policy>
    double
    sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr, const double &mu,
        const mutation_model &mmodel, const recombination_policy &rec_pol,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected, const double f,
        const mutation_removal_policy mp)
    {
        // run changing N version with N_next == N_curr
        return sample_diploid(r, gametes, diploids, mutations, mcounts, N_curr,
                              N_curr, mu, mmodel, rec_pol, ff, neutral,
                              selected, f, mp);
    }

    // single deme, N changing
    template <typename gamete_type, typename gamete_cont_type_allocator,
              typename mutation_type, typename mutation_cont_type_allocator,
              typename diploid_geno_t, typename diploid_vector_type_allocator,
              typename diploid_fitness_function, typename mutation_model,
              typename recombination_policy,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              typename mutation_removal_policy>
    double
    sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<diploid_geno_t, diploid_vector_type_allocator>
            &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr,
        const uint_t &N_next, const double &mu, const mutation_model &mmodel,
        const recombination_policy &rec_pol,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected, const double f,
        const mutation_removal_policy mp)
    {
        /*
          The main part of fwdpp does not throw exceptions.
          Rather, testing is performed via C's assert macro.
          This macro should be disabled in "production" builds via
          -DNEBUG as is standard practice.  It is the developer's
          responsibility to properly set up a build system to distinguish
          'debug' from 'production' builds.

          More complex debugging blocks will be wrapped in #ifndef
          NDEBUG/#endif
          blocks as needed.

          Compiling in a 'debug' mode slows simulations down several-fold.
        */

        // test preconditions in debugging mode
#ifndef NDEBUG
        if (mcounts.size() != mutations.size())
            {
                throw std::runtime_error(
                    "FWDPP DEBUG: mutation container size must equal "
                    "mutation count container size");
            }
        if (N_curr != diploids.size())
            {
                throw std::runtime_error(
                    "FWDPP DEBUG: N_curr != diploids.size()");
            }
#endif

        mut_rec_intermediates intermediates;
        /*
          The mutation and gamete containers contain both extinct and extant
          objects.
          The former are useful b/c the represent already-allocated memory.
          The library
          uses these extinct objects to 'recycle' them into new objects.  The
          function calls
          below create FIFO queues of where extinct objects are.  These queues
          are passed to
          mutation and recombination functions and used to decide if recyling
          is possible or
          if a new object needs to be 'emplace-back'-ed into a container.

          The type of the FIFO queue is abstracted with the name
          fwdpp::fwdpp_internal::recycling_bin_t,
          which is a C++11 template alias.

          The details of recycling are implemented in
          fwdpp/internal/recycling.hpp
        */
        intermediates.mutation_recycling_bin
            = fwdpp_internal::make_mut_queue(mcounts);
        intermediates.gamete_recycling_bin
            = fwdpp_internal::make_gamete_queue(gametes);

        // Calculate fitness for each diploid:

        // create a vector to store fitnesses:
        std::vector<double> fitnesses(diploids.size());
        double wbar = 0.; // pop'n mean fitness
        for (uint_t i = 0; i < N_curr; ++i)
            {
                /*
                  Set the count of each gamete to 0.

                  Yes, in the case of > 1 diploid having the exact same gamete
                  index,
                  this is done multiple times.
                */
                gametes[diploids[i].first].n = gametes[diploids[i].second].n
                    = 0;
                /*
                  Assign fitness to the i-th individual.

                  ff is a "fitness function", which returns a double.  For
                  examples, see
                  fwdpp::multiplicative_diploid, which is a "standard" type of
                  fitness function
                  used in population genetics.  "Standard" types of models are
                  defined in
                  fwdpp/fitness_models.hpp.
                 */
                fitnesses[i] = ff(diploids[i], gametes, mutations);
                wbar += fitnesses[i];
            }
        wbar /= double(diploids.size());
#ifndef NDEBUG
        for (const auto &g : gametes)
            {
                if (g.n > 0)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: not all gamete counts equal zero");
                    }
            }
#endif

        /*
          This is a lookup table for rapid sampling of diploids proportional to
          their fitnesses.
          This is a unique_ptr wrapper around an object from the GNU Scientific
          Library.  A custom deleter
          is required to make this work, which is why there is no cleanup call
          down below.
        */
        fwdpp_internal::gsl_ran_discrete_t_ptr lookup(
            gsl_ran_discrete_preproc(N_curr, fitnesses.data()));
        const auto parents(diploids); // Copy the parents, which is trivally
        // fast for the vast majority of use
        // cases.

        // Change the population size
        if (diploids.size() != N_next)
            {
                diploids.resize(N_next);
            }

        // Fill in the next generation!
        for (auto &dip : diploids)
            {
                // Choose parent 1 based on fitness
                auto p1 = gsl_ran_discrete(r, lookup.get());
                // If inbred (w/probability f2), parent2 = parent1, else choose
                // again based on fitness
                auto p2 = (f == 1. || (f > 0. && gsl_rng_uniform(r) < f))
                              ? p1
                              : gsl_ran_discrete(r, lookup.get());
                generate_offspring_gametes(r, mu, dip, parents[p1],
                                           parents[p2], gametes, mutations,
                                           intermediates, mmodel, rec_pol);
            }
#ifndef NDEBUG
        for (const auto &dip : diploids)
            {
                assert(gametes[dip.first].n > 0);
                assert(gametes[dip.first].n <= 2 * N_next);
                assert(gametes[dip.second].n > 0);
                assert(gametes[dip.second].n <= 2 * N_next);
            }
#endif
        /*
          At the end of the above loop, we have a bunch of new diploids
          that are all recombined and mutated sampling of the parental
          generation.

          Our problem is that we no longer know how many times each mutation is
          present, which
          is corrected by the following call.

          Although the implementation of process_gametes is super-trivial, it
          is actually the
          most computationally-expensive part of a simulation once mutation
          rates are large.

          Further, the function is hard to optimize. Recall that gametes store
          mutations in order
          according to position.  Thus, when we go from a gamete to a position
          in mcounts, we are
          accessing the latter container out of order with respect to location
          in memory.  process_gametes
          is thus the "scatter" part of a "scatter-gather" idiom.  Modern x86
          CPU have little available
          for vectorizing such cases.  I've experimented with CPU intrinsics to
          attempt memory prefetches,
          but never saw any performance improvement, and the code got complex,
          and possibly less portable.

          The implementation is in fwdpp/internal/sample_diploid_helpers.hpp
         */
        fwdpp_internal::process_gametes(gametes, mutations, mcounts);
#ifndef NDEBUG
        for (const auto &mc : mcounts)
            {
                if (mc > 2 * N_next)
                    {
                        throw std::runtime_error("mutation size too large");
                    }
            }
#endif

        /*
          The last thing to do is handle fixations.  In many contexts, we
          neither want nor need
          to keep indexes to fixed variants in our gametes.  Such decisions are
          implemented via
          simple policies, which are in the variable 'mp'.

          The implementation is in fwdpp/internal/gamete_cleaner.hpp.

          The implementation is the "erase/remove idiom" (Effective STL, Item
          32), but with a twist
          that the function will exit early if there are no fixations present
          in the population at
          the moment.

          Example policies are fwdpp::remove_nothing and fwdpp::remove_neutral,
          both found
          in fwdpp/fwd_functional.hpp.  If mp is std::true_type, then all
          fixations (e.g., neutral
          and selected)  will be removed from all gametes.
        */
        fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts, 2 * N_next,
                                       mp);
        return wbar;
    }

    // Multi-locus API
    // single deme, N changing
    template <typename diploid_geno_t, typename gamete_type,
              typename gamete_cont_type_allocator, typename mutation_type,
              typename mutation_cont_type_allocator,
              typename diploid_vector_type_allocator,
              typename locus_vector_type_allocator,
              typename diploid_fitness_function,
              typename mutation_model_container,
              typename recombination_policy_container,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              template <typename, typename> class locus_vector_type,
              typename mutation_removal_policy>
    double
    sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<
            locus_vector_type<diploid_geno_t, locus_vector_type_allocator>,
            diploid_vector_type_allocator> &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N_curr,
        const uint_t &N_next, const double *mu,
        const mutation_model_container &mmodel,
        const recombination_policy_container &rec_policies,
        const std::vector<std::function<unsigned(void)>> &interlocus_rec,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected, const double &f,
        const mutation_removal_policy &mp)
    {
#ifndef NDEBUG
        if (mcounts.size() != mutations.size())
            {
                throw std::runtime_error(
                    "FWDPP DEBUG: mutation container size must equal "
                    "mutation count container size");
            }
        if (N_curr != diploids.size())
            {
                throw std::runtime_error(
                    "FWDPP DEBUG: N_curr != diploids.size()");
            }
#endif
        // Vector of parental fitnesses
        std::vector<double> fitnesses(N_curr);
        double wbar = 0.;
        auto mut_recycling_bin = fwdpp_internal::make_mut_queue(mcounts);
        auto gamete_recycling_bin = fwdpp_internal::make_gamete_queue(gametes);
        // Go over parents
        for (uint_t i = 0; i < N_curr; ++i)
            {
                // set parental gamete counts to 0 for each locus
                for (auto locus : diploids[i])
                    {
                        gametes[locus.first].n = gametes[locus.second].n = 0;
                    }
                // Calculate the fitness of this parent
                fitnesses[i] = ff(diploids[i], gametes, mutations);
                // increment pop. mean fitness
                wbar += fitnesses[i];
            }
        wbar /= double(diploids.size());
#ifndef NDEBUG
        /*
          If we are debugging, let's make sure that every gamete has been set
          to n = 0.
          Rationale for check:  if we are failing to update data types
          properly, then
          it is possible that the "gamete pool" contains items not carried by
          any diploids.
          If so, this assertion will fail.
        */
        for (const auto &g : gametes)
            {
                if (g.n > 0)
                    {
                        throw std::runtime_error("FWDPP DEBUG: not all gamete "
                                                 "counts were set to zero");
                    }
            }
#endif

        fwdpp_internal::gsl_ran_discrete_t_ptr lookup(
            gsl_ran_discrete_preproc(fitnesses.size(), fitnesses.data()));

        const auto parents(diploids); // Copy the parents.  Exact copy of
        // diploids--same fitnesses, etc.

        // Change the population size and reset dptr to avoid iterator
        // invalidation
        if (diploids.size() != N_next)
            {
                diploids.resize(N_next);
            }

        for (auto &dip : diploids)
            {
                auto p1 = gsl_ran_discrete(r, lookup.get());
                auto p2 = (f == 1. || (f > 0. && gsl_rng_uniform(r) < f))
                              ? p1
                              : gsl_ran_discrete(r, lookup.get());
                dip = fwdpp_internal::multilocus_rec_mut(
                    r, parents[p1], parents[p2], mut_recycling_bin,
                    gamete_recycling_bin, rec_policies, interlocus_rec,
                    ((gsl_rng_uniform(r) < 0.5) ? 1 : 0),
                    ((gsl_rng_uniform(r) < 0.5) ? 1 : 0), gametes, mutations,
                    neutral, selected, mu, mmodel);
            }
        fwdpp_internal::process_gametes(gametes, mutations, mcounts);
        fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts, 2 * N_next,
                                       mp, std::true_type());
        return wbar;
    }

    // single deme, constant N
    template <typename diploid_geno_t, typename gamete_type,
              typename gamete_cont_type_allocator, typename mutation_type,
              typename mutation_cont_type_allocator,
              typename diploid_vector_type_allocator,
              typename locus_vector_type_allocator,
              typename diploid_fitness_function,
              typename mutation_model_container,
              typename recombination_policy_container,
              template <typename, typename> class gamete_cont_type,
              template <typename, typename> class mutation_cont_type,
              template <typename, typename> class diploid_vector_type,
              template <typename, typename> class locus_vector_type,
              typename mutation_removal_policy>
    double
    sample_diploid(
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        diploid_vector_type<
            locus_vector_type<diploid_geno_t, locus_vector_type_allocator>,
            diploid_vector_type_allocator> &diploids,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator>
            &mutations,
        std::vector<uint_t> &mcounts, const uint_t &N, const double *mu,
        const mutation_model_container &mmodel,
        const recombination_policy_container &rec_policies,
        const std::vector<std::function<unsigned(void)>> &interlocus_rec,
        const diploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected, const double &f,
        const mutation_removal_policy &mp)
    {
        return sample_diploid(r, gametes, diploids, mutations, mcounts, N, N,
                              mu, mmodel, rec_policies, interlocus_rec, ff,
                              neutral, selected, f, mp);
    }
} // namespace fwdpp

#endif
