/*!
  \file mutate_recombine.hpp

  \brief Handling of recombination and mutation in one step.

  \note Introduced in fwdpp 0.5.7
*/
#ifndef FWDPP_MUTATE_RECOMBINE_HPP__
#define FWDPP_MUTATE_RECOMBINE_HPP__

#include <vector>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/mutation_internal.hpp>
#include <fwdpp/internal/rec_gamete_updater.hpp>

namespace KTfwd
{
    template <typename recombination_policy, typename gcont_t,
              typename mcont_t>
    std::vector<double>
    generate_breakpoints(const std::size_t g1, const std::size_t g2,
                         const gcont_t &gametes, const mcont_t &mutations,
                         const recombination_policy &rec_pol)
    /*! Generate vector of recombination breakpoints
     \param g1 Index of gamete 1
     \param g2 Index of gamete 2
     \param gametes Vector of gametes
     \param mutation Vector of mutations
     \param rec_pol Function to generate breakpoints

     \return std::vector<double> containing sorted breakpoints

     \note An empty return value means no breakpoints.  Otherwise,
     the breakpoints are returned and are terminated by
     std::numeric_limits<double>::max()
     */
    {
        auto nm1
            = gametes[g1].mutations.size() + gametes[g1].smutations.size();
        auto nm2
            = gametes[g2].mutations.size() + gametes[g2].smutations.size();
        if ((std::min(nm1, nm2) == 0 && std::max(nm1, nm2) == 1)
            || gametes[g1] == gametes[g2])
            {
                return {};
            }
        return rec_pol(gametes[g1], gametes[g2], mutations);
    }

    template <typename queue_type, typename mutation_model, typename gcont_t,
              typename mcont_t>
    std::vector<uint_t>
    generate_new_mutations(queue_type &recycling_bin, const gsl_rng *r,
                           const double &mu, gcont_t &gametes,
                           mcont_t &mutations, const std::size_t g,
                           const mutation_model &mmodel)
    /*!
         Return a vector of keys to new mutations.  The keys
     will be sorted according to mutation postition.

         \param recycling_bin The queue for recycling mutations
         \param r A random number generator
         \param mu The total mutation rate
		 \param gametes Vector of gametes
         \param mutations Vector of mutations
		 \param g index of gamete to mutate
         \param mmodel The mutation policy

         \return Vector of mutation keys, sorted according to position
	 */
    {
        unsigned nm = gsl_ran_poisson(r, mu);
        std::vector<uint_t> rv;
        rv.reserve(nm);
        for (unsigned i = 0; i < nm; ++i)
            {
                rv.emplace_back(fwdpp_internal::mmodel_dispatcher(
                    mmodel, gametes[g], mutations, recycling_bin));
            }
        std::sort(rv.begin(), rv.end(),
                  [&mutations](const uint_t a, const uint_t b) {
                      return mutations[a].pos < mutations[b].pos;
                  });
        return rv;
    }

    namespace fwdpp_internal
    {
        template <typename container, typename integer_type, typename mcont_t>
        typename container::iterator
        insert_new_mutation(const typename container::iterator beg,
                            const typename container::iterator end,
                            const integer_type mut_key,
                            const mcont_t &mutations, container &c)
        //Inserts mutation key into c such that sort order is maintained
        {
            auto t = std::upper_bound(
                beg, end, mutations[mut_key].pos,
                [&mutations](const double &v, const uint_t mut) noexcept {
                    return v < mutations[mut].pos;
                });
            c.insert(c.end(), beg, t);
            c.push_back(mut_key);
            return t;
        }

        template <typename gcont_t, typename container>
        void
        prep_temporary_containers(const std::size_t g1, const std::size_t g2,
                                  const gcont_t &gametes, container &neutral,
                                  container &selected)
        //Clear temporary containers and reserve memory
        {
            neutral.clear();
            selected.clear();
            neutral.reserve(std::max(gametes[g1].mutations.size(),
                                     gametes[g2].mutations.size()));
            selected.reserve(std::max(gametes[g1].smutations.size(),
                                      gametes[g2].smutations.size()));
        }
    }

    template <typename gcont_t, typename mcont_t, typename queue_type>
    uint_t
    mutate_recombine(
        const std::vector<uint_t> &new_mutations,
        const std::vector<double> &breakpoints, const std::size_t g1,
        const std::size_t g2, gcont_t &gametes, mcont_t &mutations,
        queue_type &gamete_recycling_bin,
        typename gcont_t::value_type::mutation_container &neutral,
        typename gcont_t::value_type::mutation_container &selected)
    {
        if (new_mutations.empty() && breakpoints.empty())
            {
                return g1;
            }
        else if (breakpoints.empty()) // only mutations to deal with
            {
                fwdpp_internal::prep_temporary_containers(g1, g2, gametes,
                                                          neutral, selected);
                auto nb = gametes[g1].mutations.begin(),
                     sb = gametes[g1].smutations.begin();
                const auto ne = gametes[g1].mutations.end(),
                           se = gametes[g1].smutations.end();
                for (auto &&m : new_mutations)
                    {
                        if (mutations[m].neutral)
                            {
                                nb = fwdpp_internal::insert_new_mutation(
                                    nb, ne, m, mutations, neutral);
                            }
                        else
                            {
                                sb = fwdpp_internal::insert_new_mutation(
                                    sb, se, m, mutations, selected);
                            }
                    }
                neutral.insert(neutral.end(), nb, ne);
                selected.insert(selected.end(), sb, se);

                return fwdpp_internal::recycle_gamete(
                    gametes, gamete_recycling_bin, neutral, selected);
            }
        // If we get here, there are mutations and
        // recombinations to handle
        fwdpp_internal::prep_temporary_containers(g1, g2, gametes, neutral,
                                                  selected);

        auto itr = gametes[g1].mutations.cbegin();
        auto jtr = gametes[g2].mutations.cbegin();
        auto itr_s = gametes[g1].smutations.cbegin();
        auto jtr_s = gametes[g2].smutations.cbegin();
        auto itr_e = gametes[g1].mutations.cend();
        auto itr_s_e = gametes[g1].smutations.cend();
        auto jtr_e = gametes[g2].mutations.cend();
        auto jtr_s_e = gametes[g2].smutations.cend();

        auto next_mutation = new_mutations.cbegin();
        for (auto i = breakpoints.cbegin(); i != breakpoints.cend();)
            {
                if (next_mutation != new_mutations.cend()
                    && mutations[*next_mutation].pos < *i)
                    {
                        const auto mut = &mutations[*next_mutation];
                        itr = fwdpp_internal::rec_gam_updater(
                            itr, itr_e, mutations, neutral, mut->pos);
                        itr_s = fwdpp_internal::rec_gam_updater(
                            itr_s, itr_s_e, mutations, selected, mut->pos);
                        jtr = fwdpp_internal::rec_update_itr(
                            jtr, jtr_e, mutations, mut->pos);
                        jtr_s = fwdpp_internal::rec_update_itr(
                            jtr_s, jtr_s_e, mutations, mut->pos);
                        if (mut->neutral)
                            {
                                neutral.push_back(*next_mutation);
                            }
                        else
                            {
                                selected.push_back(*next_mutation);
                            }
                        ++next_mutation;
                    }
                else
                    {
                        itr = fwdpp_internal::rec_gam_updater(
                            itr, itr_e, mutations, neutral, *i);
                        itr_s = fwdpp_internal::rec_gam_updater(
                            itr_s, itr_s_e, mutations, selected, *i);
                        jtr = fwdpp_internal::rec_update_itr(jtr, jtr_e,
                                                             mutations, *i);
                        jtr_s = fwdpp_internal::rec_update_itr(jtr_s, jtr_s_e,
                                                               mutations, *i);
                        std::swap(itr, jtr);
                        std::swap(itr_s, jtr_s);
                        std::swap(itr_e, jtr_e);
                        std::swap(itr_s_e, jtr_s_e);
                        ++i;
                    }
            }
        assert(next_mutation == new_mutations.cend());
        return fwdpp_internal::recycle_gamete(gametes, gamete_recycling_bin,
                                              neutral, selected);
    }
}

#endif
