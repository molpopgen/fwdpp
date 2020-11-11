/*!
  \file mutate_recombine.hpp

  \brief Handling of recombination and mutation in one step.

  \note Introduced in fwdpp 0.5.7
*/
#ifndef FWDPP_MUTATE_RECOMBINE_HPP__
#define FWDPP_MUTATE_RECOMBINE_HPP__

#include <vector>
#include <algorithm>
#include <tuple>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/debug.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/mutation_internal.hpp>
#include <fwdpp/internal/rec_gamete_updater.hpp>
#include <fwdpp/internal/recombination_common.hpp>

namespace fwdpp
{
    template <typename recombination_policy, typename DiploidType,
              typename haploid_genome_t, typename MutationContainerType>
    inline typename std::result_of<recombination_policy()>::type
    dispatch_recombination_policy(const recombination_policy &rec_pol, DiploidType &&,
                                  haploid_genome_t &&, haploid_genome_t &&,
                                  MutationContainerType &&)
    {
        return rec_pol();
    }

    template <typename recombination_policy, typename DiploidType,
              typename haploid_genome_t, typename MutationContainerType>
    inline typename std::result_of<recombination_policy(
        haploid_genome_t &&, haploid_genome_t &&, MutationContainerType &&)>::type
    dispatch_recombination_policy(const recombination_policy &rec_pol, DiploidType &&,
                                  haploid_genome_t &&g1, haploid_genome_t &&g2,
                                  MutationContainerType &&mutations)
    {
        return rec_pol(std::forward<haploid_genome_t>(g1),
                       std::forward<haploid_genome_t>(g2),
                       std::forward<MutationContainerType>(mutations));
    }

    template <typename recombination_policy, typename DiploidType,
              typename haploid_genome_t, typename MutationContainerType>
    inline
        typename std::result_of<recombination_policy(DiploidType &&, haploid_genome_t &&,
                                                     haploid_genome_t &&,
                                                     MutationContainerType &&)>::type
        dispatch_recombination_policy(const recombination_policy &rec_pol,
                                      DiploidType &&diploid, haploid_genome_t &&g1,
                                      haploid_genome_t &&g2,
                                      MutationContainerType &&mutations)
    {
        return rec_pol(std::forward<DiploidType>(diploid),
                       std::forward<haploid_genome_t>(g1),
                       std::forward<haploid_genome_t>(g2),
                       std::forward<MutationContainerType>(mutations));
    }

    template <typename recombination_policy, typename DiploidType,
              typename GenomeContainerType, typename MutationContainerType>
    std::vector<double>
    generate_breakpoints(const DiploidType &diploid, const std::size_t g1,
                         const std::size_t g2,
                         const GenomeContainerType &haploid_genomes,
                         const MutationContainerType &mutations,
                         const recombination_policy &rec_pol)
    /// Generate vector of recombination breakpoints
    ///
    /// \param diploid A single-locus diploid.
    /// \param g1 Index of haploid_genome 1
    /// \param g2 Index of haploid_genome 2
    /// \param haploid_genomes Vector of haploid_genomes
    /// \param mutations Vector of mutations
    /// \param rec_pol Function to generate breakpoints
    ///
    /// \return std::vector<double> containing sorted breakpoints
    ///
    /// \note An empty return value means no breakpoints.  Otherwise,
    /// the breakpoints are returned and are terminated by
    /// std::numeric_limits<double>::max()
    {
        //TODO: decide if we wish to re-enable the code below.
        //auto nm1
        //    = haploid_genomes[g1].mutations.size() + haploid_genomes[g1].smutations.size();
        //auto nm2
        //    = haploid_genomes[g2].mutations.size() + haploid_genomes[g2].smutations.size();
        //if ((std::min(nm1, nm2) == 0 && std::max(nm1, nm2) == 1)
        //    || haploid_genomes[g1] == haploid_genomes[g2])
        //    {
        //        return {};
        //    }
        return dispatch_recombination_policy(
            std::cref(rec_pol), std::cref(diploid), std::cref(haploid_genomes[g1]),
            std::cref(haploid_genomes[g2]), std::cref(mutations));
    }

    template <typename mutation_model, typename DiploidType, typename GenomeContainerType,
              typename MutationContainerType>
    std::vector<uint_t>
    generate_new_mutations(flagged_mutation_queue &recycling_bin, const gsl_rng *r,
                           const double &mu, const DiploidType &dip,
                           GenomeContainerType &haploid_genomes,
                           MutationContainerType &mutations, const std::size_t g,
                           const mutation_model &mmodel)
    ///
    /// Return a vector of keys to new mutations.  The keys
    /// will be sorted according to mutation postition.
    ///
    /// \param recycling_bin The queue for recycling mutations
    /// \param r A random number generator
    /// \param mu The total mutation rate
    /// \param dip A single-locus diploid
    /// \param haploid_genomes Vector of haploid_genomes
    /// \param mutations Vector of mutations
    /// \param g index of haploid_genome to mutate
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
                    mmodel, dip, haploid_genomes[g], mutations, recycling_bin));
            }
        std::sort(rv.begin(), rv.end(), [&mutations](const uint_t a, const uint_t b) {
            return mutations[a].pos < mutations[b].pos;
        });
        return rv;
    }

    namespace fwdpp_internal
    {
        template <typename container, typename integer_type,
                  typename MutationContainerType>
        typename container::iterator
        insert_new_mutation(const typename container::iterator beg,
                            const typename container::iterator end,
                            const integer_type mut_key,
                            const MutationContainerType &mutations, container &c)
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

        template <typename GenomeContainerType, typename container>
        void
        prep_temporary_containers(const std::size_t g1, const std::size_t g2,
                                  const GenomeContainerType &haploid_genomes,
                                  container &neutral, container &selected)
        // Clear temporary containers and reserve memory
        {
            neutral.clear();
            selected.clear();
            neutral.reserve(std::max(haploid_genomes[g1].mutations.size(),
                                     haploid_genomes[g2].mutations.size()));
            selected.reserve(std::max(haploid_genomes[g1].smutations.size(),
                                      haploid_genomes[g2].smutations.size()));
        }
    } // namespace fwdpp_internal

    template <typename GenomeContainerType, typename MutationContainerType,
              typename queue_type, typename new_mutations_type,
              typename breakpoints_type>
    uint_t
    mutate_recombine(
        const new_mutations_type &new_mutations, const breakpoints_type &breakpoints,
        const std::size_t g1, const std::size_t g2, GenomeContainerType &haploid_genomes,
        MutationContainerType &mutations, queue_type &haploid_genome_recycling_bin,
        typename GenomeContainerType::value_type::mutation_container &neutral,
        typename GenomeContainerType::value_type::mutation_container &selected)
    ///
    /// Update apply new mutations and recombination events to
    /// an offspring's haploid_genome.
    ///
    /// \param new_mutations A range referring to new mutation keys.  Must be traversable via non-member begin/end
    /// functions using argument-dependent lookup (ADL).
    /// \param breakpoints A range of recombination breakpoints traversable via non-member begin/end functions.
    /// \param g1 Parental haploid_genome 1
    /// \param g2 Parental haploid_genome 2
    /// \param haploid_genomes The vector of haploid_genomes in the population
    /// \param mutations The vector of mutations in the population
    /// \param haploid_genome_recycling_bin FIFO queue for haploid_genome recycling
    /// \param neutral Temporary container for updating neutral mutations
    /// \param selected Temporary container for updatng selected positions
    ///
    /// \return The index of the new offspring haploid_genome in \a haploid_genomes.
    ///
    /// TODO: unit test this note.  I simply cannot believe it! :)
    /// \note For efficiency, it is helpful if \a new_mutations is sorted
    /// by mutation position.  fwdpp::generate_new_mutations exists to help in
    /// that
    /// regard. Many of the evolve functions used in this library and other
    /// packages by the author will use fwdpp::generate_breakpoints to
    /// generate \a breakpoints.  That is not, however, required.
    ///
    /// \todo Need a unit test on what happens if \a mutation_keys is not a sorted range
    ///
    /// \version
    /// This function was added in fwdpp 0.5.7.
    {
        if (begin(new_mutations) == end(new_mutations)
            && begin(breakpoints) == end(breakpoints))
            {
                return g1;
            }
        else if (begin(breakpoints) == end(breakpoints)) // only mutations to deal with
            {
                fwdpp_internal::prep_temporary_containers(g1, g2, haploid_genomes,
                                                          neutral, selected);
                auto nb = haploid_genomes[g1].mutations.begin(),
                     sb = haploid_genomes[g1].smutations.begin();
                const auto ne = haploid_genomes[g1].mutations.end(),
                           se = haploid_genomes[g1].smutations.end();
                for (auto m = begin(new_mutations); m < end(new_mutations); ++m)
                    {
                        if (mutations[*m].neutral)
                            {
                                nb = fwdpp_internal::insert_new_mutation(
                                    nb, ne, *m, mutations, neutral);
                            }
                        else
                            {
                                sb = fwdpp_internal::insert_new_mutation(
                                    sb, se, *m, mutations, selected);
                            }
                    }
                neutral.insert(neutral.end(), nb, ne);
                selected.insert(selected.end(), sb, se);

#ifndef NDEBUG
                std::size_t new_neutral = 0, new_selected = 0;
                for (auto m : new_mutations)
                    {
                        if (mutations[m].neutral)
                            {
                                ++new_neutral;
                            }
                        if (!mutations[m].neutral)
                            {
                                ++new_selected;
                            }
                    }
                if (neutral.size() != haploid_genomes[g1].mutations.size() + new_neutral)
                    {
                        throw std::runtime_error("FWDPP DEBUG: failure to add "
                                                 "all new neutral mutations");
                    }
                if (selected.size()
                    != haploid_genomes[g1].smutations.size() + new_selected)
                    {
                        throw std::runtime_error("FWDPP DEBUG: failure to add "
                                                 "all new selected mutations");
                    }
#endif
                return recycle_haploid_genome(
                    haploid_genomes, haploid_genome_recycling_bin, neutral, selected);
            }
        if (breakpoints.size() == 1)
            {
                throw std::runtime_error(
                    "invalid number of breakpoints. likely sentinel error");
            }
        // If we get here, there are mutations and
        // recombinations to handle
        fwdpp_internal::prep_temporary_containers(g1, g2, haploid_genomes, neutral,
                                                  selected);

        auto itr = haploid_genomes[g1].mutations.cbegin();
        auto jtr = haploid_genomes[g2].mutations.cbegin();
        auto itr_s = haploid_genomes[g1].smutations.cbegin();
        auto jtr_s = haploid_genomes[g2].smutations.cbegin();
        auto itr_e = haploid_genomes[g1].mutations.cend();
        auto itr_s_e = haploid_genomes[g1].smutations.cend();
        auto jtr_e = haploid_genomes[g2].mutations.cend();
        auto jtr_s_e = haploid_genomes[g2].smutations.cend();

        auto next_mutation = begin(new_mutations);
        for (auto i = begin(breakpoints); i != end(breakpoints);)
            {
                if (next_mutation != end(new_mutations)
                    && mutations[*next_mutation].pos < *i)
                    {
                        const auto mut = &mutations[*next_mutation];
                        itr = fwdpp_internal::rec_gam_updater(itr, itr_e, mutations,
                                                              neutral, mut->pos);
                        itr_s = fwdpp_internal::rec_gam_updater(
                            itr_s, itr_s_e, mutations, selected, mut->pos);
                        jtr = fwdpp_internal::rec_update_itr(jtr, jtr_e, mutations,
                                                             mut->pos);
                        jtr_s = fwdpp_internal::rec_update_itr(jtr_s, jtr_s_e, mutations,
                                                               mut->pos);
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
                        itr = fwdpp_internal::rec_gam_updater(itr, itr_e, mutations,
                                                              neutral, *i);
                        itr_s = fwdpp_internal::rec_gam_updater(itr_s, itr_s_e,
                                                                mutations, selected, *i);
                        jtr = fwdpp_internal::rec_update_itr(jtr, jtr_e, mutations, *i);
                        jtr_s = fwdpp_internal::rec_update_itr(jtr_s, jtr_s_e, mutations,
                                                               *i);
                        std::swap(itr, jtr);
                        std::swap(itr_s, jtr_s);
                        std::swap(itr_e, jtr_e);
                        std::swap(itr_s_e, jtr_s_e);
                        ++i;
                    }
            }
#ifndef NDEBUG
        if (next_mutation != end(new_mutations))
            {
                throw std::runtime_error("FWDPP DEBUG: fatal error during "
                                         "mutation/recombination");
            }
#endif
        return recycle_haploid_genome(haploid_genomes, haploid_genome_recycling_bin,
                                      neutral, selected);
    }

    template <typename DiploidType, typename GenomeContainerType,
              typename MutationContainerType, typename recmodel, typename mutmodel>
    std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>
    mutate_recombine_update(
        const gsl_rng *r, GenomeContainerType &haploid_genomes,
        MutationContainerType &mutations,
        std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>
            parental_haploid_genomes,
        const recmodel &rec_pol, const mutmodel &mmodel, const double mu,
        flagged_haploid_genome_queue &haploid_genome_recycling_bin,
        flagged_mutation_queue &mutation_recycling_bin, DiploidType &dip,
        typename GenomeContainerType::value_type::mutation_container &neutral,
        typename GenomeContainerType::value_type::mutation_container &selected)
    ///
    /// "Convenience" function for generating offspring haploid_genomes.
    ///
    /// This function calls fwdpp::generate_breakpoints,
    /// fwdpp::generate_new_mutations, and fwdpp::mutate_recombine,
    /// resulting in two offspring haploid_genomes.
    ///
    /// \param r A gsl_rng *
    /// \param haploid_genomes Vector of haploid_genomes in population
    /// \param mutations Vector of mutations in population
    /// \param parental_haploid_genomes Tuple of haploid_genome keys for each parent
    /// \param rec_pol Policy to generate recombination breakpoints
    /// \param mmodel Policy to generate new mutations
    /// \param mu Total mutation rate (per haploid_genome).
    /// \param haploid_genome_recycling_bin FIFO queue for haploid_genome recycling
    /// \param mutation_recycling_bin FIFO queue for mutation recycling
    /// \param dip The offspring
    /// \param neutral Temporary container for updating neutral mutations
    /// \param selected Temporary container for updating selected mutations
    ///
    /// \return Number of recombination breakpoints and mutations for each
    /// haploid_genome.
    ///
    /// \note
    /// \a parental_haploid_genomes should contain parent one/haploid_genome one,
    /// parent one/haploid_genome two, parent two/haploid_genome one,
    /// and parent two/haploid_genome two, in that order.
    ///
    /// \version
    /// Added in fwdpp 0.5.7.
    {
        auto p1g1 = std::get<0>(parental_haploid_genomes);
        auto p1g2 = std::get<1>(parental_haploid_genomes);
        auto p2g1 = std::get<2>(parental_haploid_genomes);
        auto p2g2 = std::get<3>(parental_haploid_genomes);
        // Now, we generate the crossover breakpoints for
        // both parents,as well as the new mutations that we'll place
        // onto each haploid_genome.  The specific order of operations below
        // is done to ensure the exact same output as fwdpp 0.5.6 and
        // earlier.
        // The breakpoints are of type std::vector<double>, and
        // the new_mutations are std::vector<fwdpp::uint_t>, with
        // the integers representing the locations of the new mutations
        // in "mutations".

        auto breakpoints
            = generate_breakpoints(dip, p1g1, p1g2, haploid_genomes, mutations, rec_pol);
        auto breakpoints2
            = generate_breakpoints(dip, p2g1, p2g2, haploid_genomes, mutations, rec_pol);
        auto new_mutations
            = generate_new_mutations(mutation_recycling_bin, r, mu, dip, haploid_genomes,
                                     mutations, p1g1, mmodel);
        auto new_mutations2
            = generate_new_mutations(mutation_recycling_bin, r, mu, dip, haploid_genomes,
                                     mutations, p2g1, mmodel);

        // Pass the breakpoints and new mutation keys on to
        // fwdpp::mutate_recombine (defined in
        // fwdpp/mutate_recombine.hpp),
        // which splices together the offspring haploid_genome and returns its
        // location in haploid_genomes.  The location of the offspring haploid_genome
        // is
        // either reycled from an extinct haploid_genome or it is the location
        // of a
        // new haploid_genome emplace_back'd onto the end.
        dip.first = mutate_recombine(new_mutations, breakpoints, p1g1, p1g2,
                                     haploid_genomes, mutations,
                                     haploid_genome_recycling_bin, neutral, selected);
        debug::haploid_genome_is_sorted(haploid_genomes[dip.first], mutations);
        dip.second = mutate_recombine(new_mutations2, breakpoints2, p2g1, p2g2,
                                      haploid_genomes, mutations,
                                      haploid_genome_recycling_bin, neutral, selected);
        debug::haploid_genome_is_sorted(haploid_genomes[dip.second], mutations);
        haploid_genomes[dip.first].n++;
        haploid_genomes[dip.second].n++;

        debug::haploid_genome_is_extant(haploid_genomes[dip.first]);
        debug::haploid_genome_is_extant(haploid_genomes[dip.second]);

        auto nrec = (!breakpoints.empty()) ? breakpoints.size() - 1 : 0;
        auto nrec2 = (!breakpoints2.empty()) ? breakpoints2.size() - 1 : 0;
        return std::make_tuple(nrec, nrec2, new_mutations.size(), new_mutations2.size());
    }
} // namespace fwdpp

#endif
