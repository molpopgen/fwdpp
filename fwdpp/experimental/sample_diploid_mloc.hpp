/*!
  \file fwdpp/experimental/sample_diploid_mloc.hpp
  \brief Testing ground for more flexible API to evolve populations --
  multilocus version.
*/
#ifndef __FWDPP_EXPERIMENTAL_SAMPLE_DIPLOID_MLOC_HPP__
#define __FWDPP_EXPERIMENTAL_SAMPLE_DIPLOID_MLOC_HPP__

#include <fwdpp/experimental/dispatch.hpp>

namespace KTfwd
{
    namespace experimental
    {
        /*!
          \brief Abstraction of the standard Wright-Fisher sampling process for
          a multi-locus/region simulation
        */
        struct standardWFrules_mloc
        {
            double wbar;
            std::vector<double> fitnesses;

            fwdpp_internal::gsl_ran_discrete_t_ptr lookup;
            //! \brief Constructor
            standardWFrules_mloc()
                : wbar(0.), fitnesses(std::vector<double>()),
                  lookup(fwdpp_internal::gsl_ran_discrete_t_ptr(nullptr))
            {
            }

            //! \brief The "fitness manager"
            template <typename dipcont_t, typename gcont_t, typename mcont_t,
                      typename fitness_func>
            void
            w(const dipcont_t &diploids, gcont_t &gametes,
              const mcont_t &mutations, const fitness_func &ff)
            {
                using diploid_geno_t = typename dipcont_t::value_type;
                unsigned N_curr = diploids.size();
                if (fitnesses.size() < N_curr)
                    fitnesses.resize(N_curr);
                wbar = 0.;

                for (unsigned i = 0; i < N_curr; ++i)
                    {
                        for (auto region : diploids[i])
                            {
                                gametes[region.first].n
                                    = gametes[region.second].n = 0;
                            }
                        // No need to dispatch--this type of simulation is
                        // always
                        // regarded as a non-standard diploid type.
                        fitnesses[i] = ff(diploids[i], gametes, mutations);
                    }

                wbar /= double(diploids.size());

                /*!
                  Black magic alert:
                  fwdpp_internal::gsl_ran_discrete_t_ptr contains a
                  std::unique_ptr wrapping the GSL pointer.
                  This type has its own deleter, which is convenient, because
                  operator= for unique_ptrs automagically calls the deleter
                  before assignment!
                  Details:
                  http://www.cplusplus.com/reference/memory/unique_ptr/operator=

                  This only works b/c the rhs of the expression below may be
                  treated as an rvalue reference.
                */
                lookup = fwdpp_internal::gsl_ran_discrete_t_ptr(
                    gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
            }

            //! \brief Pick parent one
            inline size_t
            pick1(const gsl_rng *r) const
            {
                return gsl_ran_discrete(r, lookup.get());
            }

            //! \brief Pick parent 2.  Parent 1's data are passed along for
            //! models where that is relevant
            template <typename diploid_t, typename gcont_t, typename mcont_t>
            inline size_t
            pick2(const gsl_rng *r, const size_t &p1, const double &f,
                  const diploid_t &, const gcont_t &, const mcont_t &) const
            {
                return ((f == 1.) || (f > 0. && gsl_rng_uniform(r) < f))
                           ? p1
                           : gsl_ran_discrete(r, lookup.get());
            }

            //! \brief Update some property of the offspring based on
            //! properties of the parents
            template <typename diploid_t, typename gcont_t, typename mcont_t>
            void
            update(const gsl_rng *, diploid_t &, const diploid_t &,
                   const diploid_t &, const gcont_t &, const mcont_t &) const
            {
            }
        };

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
                  typename rules_type = standardWFrules_mloc,
                  typename mutation_removal_policy = std::true_type,
                  typename gamete_insertion_policy = emplace_back>
        double
        sample_diploid(
            const gsl_rng *r,
            gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
            diploid_vector_type<locus_vector_type<diploid_geno_t,
                                                  locus_vector_type_allocator>,
                                diploid_vector_type_allocator> &diploids,
            mutation_cont_type<mutation_type, mutation_cont_type_allocator>
                &mutations,
            std::vector<uint_t> &mcounts, const uint_t &N_curr,
            const uint_t &N_next, const double *mu,
            const mutation_model_container &mmodel,
            const recombination_policy_container &rec_policies,
			const std::vector<std::function<unsigned(void)>> & interlocus_rec,
            const diploid_fitness_function &ff,
            typename gamete_type::mutation_container &neutral,
            typename gamete_type::mutation_container &selected,
            const double &f = 0, rules_type &&rules = rules_type(),
            const mutation_removal_policy &mp = mutation_removal_policy(),
            const gamete_insertion_policy &gpolicy_mut
            = gamete_insertion_policy())
        /*
          Changing N over time
        */
        {
            auto gamete_recycling_bin
                = fwdpp_internal::make_gamete_queue(gametes);
            auto mutation_recycling_bin
                = fwdpp_internal::make_mut_queue(mcounts);

            dispatch_w(rules, diploids, gametes, mutations, ff);

            const auto parents(diploids); // copy the parents

            // Change the population size
            if (diploids.size() != N_next)
                {
                    diploids.resize(N_next);
                }
            for (auto &dip : diploids)
                {
                    // Pick parent 1
                    size_t p1 = rules.pick1(r);
                    // Pick parent 2
                    size_t p2 = rules.pick2(r, p1, f, parents[p1], gametes,
                                            mutations);
                    assert(p1 < parents.size());
                    assert(p2 < parents.size());

                    dip = fwdpp_internal::multilocus_rec_mut(
                        r, parents[p1], parents[p2], mutation_recycling_bin,
                        gamete_recycling_bin, rec_policies, interlocus_rec,
                        ((gsl_rng_uniform(r) < 0.5) ? 1 : 0),
                        ((gsl_rng_uniform(r) < 0.5) ? 1 : 0), gametes,
                        mutations, neutral, selected, mu, mmodel, gpolicy_mut);
                    dispatch_update(rules, r, dip, parents[p1], parents[p2],
                                    gametes, mutations, ff);
                }
            fwdpp_internal::process_gametes(gametes, mutations, mcounts);
            assert(popdata_sane_multilocus(diploids, gametes, mutations, mcounts));
            fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts,
                                           2 * N_next, mp, std::true_type());
            assert(check_sum(gametes, 2 * diploids[0].size() * N_next));
            return rules.wbar;
        }

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
                  typename rules_type = standardWFrules_mloc,
                  typename mutation_removal_policy = std::true_type,
                  typename gamete_insertion_policy = emplace_back>
        double
        sample_diploid(
            const gsl_rng *r,
            gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
            diploid_vector_type<locus_vector_type<diploid_geno_t,
                                                  locus_vector_type_allocator>,
                                diploid_vector_type_allocator> &diploids,
            mutation_cont_type<mutation_type, mutation_cont_type_allocator>
                &mutations,
            std::vector<uint_t> &mcounts, const uint_t &N, const double *mu,
            const mutation_model_container &mmodel,
            const recombination_policy_container &rec_policies,
            const std::vector<std::function<unsigned(void)>> &interlocus_rec,
            const diploid_fitness_function &ff,
            typename gamete_type::mutation_container &neutral,
            typename gamete_type::mutation_container &selected,
            const double &f = 0, rules_type &&rules = rules_type(),
            const mutation_removal_policy &mp = mutation_removal_policy(),
            const gamete_insertion_policy &gpolicy_mut
            = gamete_insertion_policy())
        /*
          N constant
        */
        {
            return experimental::sample_diploid(
                r, gametes, diploids, mutations, mcounts, N, N, mu, mmodel,
                rec_policies, interlocus_rec, ff, neutral, selected, f, rules,
                mp, gpolicy_mut);
        }
    }
}

#endif
