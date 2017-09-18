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
#include <tuple>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/type_traits.hpp>
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
    /// Generate vector of recombination breakpoints
    ///
    /// \param g1 Index of gamete 1
    /// \param g2 Index of gamete 2
    /// \param gametes Vector of gametes
    /// \param mutation Vector of mutations
    /// \param rec_pol Function to generate breakpoints
    ///
    /// \return std::vector<double> containing sorted breakpoints
    ///
    /// \note An empty return value means no breakpoints.  Otherwise,
    /// the breakpoints are returned and are terminated by
    /// std::numeric_limits<double>::max()
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
    ///
    /// Return a vector of keys to new mutations.  The keys
    /// will be sorted according to mutation postition.
    ///
    /// \param recycling_bin The queue for recycling mutations
    /// \param r A random number generator
    /// \param mu The total mutation rate
    /// \param gametes Vector of gametes
    /// \param mutations Vector of mutations
    /// \param g index of gamete to mutate
    /// \param mmodel The mutation policy
    ///
    /// \return Vector of mutation keys, sorted according to position
    ///
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
        // Inserts mutation key into c such that sort order is maintained
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
        // Clear temporary containers and reserve memory
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
    ///
    /// Update apply new mutations and recombination events to 
    /// an offspring's gamete.
    ///
    /// \param new_mutations Keys to new mutations
    /// \param breakpoints Recombination breakpoints
    /// \param g1 Parental gamete 1
    /// \param g2 Parental gamete 2
    /// \param gametes The vector of gametes in the population
    /// \param mutation The vector of mutations in the population
    /// \param gamete_recycling_bin FIFO queue for gamete recycling
    /// \param neutral Temporary container for updating neutral mutations
    /// \param selected Temporary container for updatng selected positions
    ///
    /// \return The index of the new offspring gamete in \a gametes.
    ///
    /// \note For efficiency, it is helpful if \a new_mutations is sorted
    /// by mutation position.  KTfwd::generate_new_mutations exists to help in that
    /// regard. Many of the evolve functions used in this library and other
    /// packages by the author will use KTfwd::generate_breakpoints to 
    /// generate \a breakpoints.  That is not, however, required.
    /// 
    /// \version
    /// This function was added in fwdpp 0.5.7.
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

    template <typename diploid_t, typename gcont_t, typename mcont_t,
              typename recmodel, typename mutmodel>
    std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>
    mutate_recombine_update(
        const gsl_rng *r, gcont_t &gametes, mcont_t &mutations,
        std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>
            parental_gametes,
        const recmodel &rec_pol, const mutmodel &mmodel, const double mu,
        typename traits::recycling_bin_t<gcont_t> &gamete_recycling_bin,
        typename traits::recycling_bin_t<mcont_t> &mutation_recycling_bin,
        diploid_t &dip,
        typename gcont_t::value_type::mutation_container &neutral,
        typename gcont_t::value_type::mutation_container &selected)
    ///
    /// "Convenience" function for generating offspring gametes.  
    /// 
    /// This function calls KTfwd::generate_breakpoints, 
    /// KTfwd::generate_new_mutations, and KTfwd::mutate_recombine,
    /// resulting in two offspring gametes.
    ///
    /// \param r A gsl_rng *
    /// \param gametes Vector of gametes in population
    /// \param mutations Vector of mutations in population
    /// \param parental_gametes Tuple of gamete keys for each parent
    /// \param rec_pol Policy to generate recombination breakpoints
    /// \param mmodel Policy to generate new mutations
    /// \param gamete_recycling_bin FIFO queue for gamete recycling
    /// \param mutation_recycling_bin FIFO queue for mutation recycling
    /// \param dip The offspring
    /// \param neutral Temporary container for updating neutral mutations
    /// \param selected Temporary container for updating selected mutations
    ///
    /// \return Number of recombination breakpoints and mutations for each
    /// gamete.
    ///
    /// \note
    /// \a parental_gametes should contain parent one/gamete one, 
    /// parent one/gamete two, parent two/gamete one, 
    /// and parent two/gamete two, in that order.
    /// 
    /// \version
    /// Added in fwdpp 0.5.7.
    {
        auto p1g1 = std::get<0>(parental_gametes);
        auto p1g2 = std::get<1>(parental_gametes);
        auto p2g1 = std::get<2>(parental_gametes);
        auto p2g2 = std::get<3>(parental_gametes);
        // Now, we generate the crossover breakpoints for
        // both parents,as well as the new mutations that we'll place
        // onto each gamete.  The specific order of operations below
        // is done to ensure the exact same output as fwdpp 0.5.6 and
        // earlier.
        // The breakpoints are of type std::vector<double>, and
        // the new_mutations are std::vector<KTfwd::uint_t>, with
        // the integers representing the locations of the new mutations
        // in "mutations".

        auto breakpoints
            = generate_breakpoints(p1g1, p1g2, gametes, mutations, rec_pol);
        auto breakpoints2
            = generate_breakpoints(p2g1, p2g2, gametes, mutations, rec_pol);
        auto new_mutations = generate_new_mutations(
            mutation_recycling_bin, r, mu, gametes, mutations, p1g1, mmodel);
        auto new_mutations2 = generate_new_mutations(
            mutation_recycling_bin, r, mu, gametes, mutations, p2g1, mmodel);

        // Pass the breakpoints and new mutation keys on to
        // KTfwd::mutate_recombine (defined in
        // fwdpp/mutate_recombine.hpp),
        // which splices together the offspring gamete and returns its
        // location in gametes.  The location of the offspring gamete
        // is
        // either reycled from an extinct gamete or it is the location
        // of a
        // new gamete emplace_back'd onto the end.
        dip.first = mutate_recombine(new_mutations, breakpoints, p1g1, p1g2,
                                     gametes, mutations, gamete_recycling_bin,
                                     neutral, selected);
        dip.second = mutate_recombine(new_mutations2, breakpoints2, p2g1, p2g2,
                                      gametes, mutations, gamete_recycling_bin,
                                      neutral, selected);
        gametes[dip.first].n++;
        gametes[dip.second].n++;

        assert(gametes[dip.first].n);
        assert(gametes[dip.second].n);

        auto nrec = (!breakpoints.empty()) ? breakpoints.size() - 1 : 0;
        auto nrec2 = (!breakpoints2.empty()) ? breakpoints2.size() - 1 : 0;
        return std::make_tuple(nrec, nrec2, new_mutations.size(),
                               new_mutations2.size());
    }
}

#endif
