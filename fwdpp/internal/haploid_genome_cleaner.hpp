#ifndef __FWDPP_INTERNAL_HAPLOID_GENOME_CLEANER_HPP__
#define __FWDPP_INTERNAL_HAPLOID_GENOME_CLEANER_HPP__

#include <vector>
#include <algorithm>
#include <limits>
#include <utility>
#include <type_traits>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>

/*!
  \file haploid_genome_cleaner.hpp

  This file handles the "pruning" of fixations from haploid_genomes.

  Whether or not any pruning happens depends on policies.

  Version 0.4.7 of the library made a major change resulting in improved
  run-time performance, especially for low mutation rates.

  Consider the case of neutral mutations occurring at rate mu (per site,
  per generation).  The rate of fixation is also mu.  Thus, when mu is low,
  the rate of fixation is low, and the expected time in generations between
  fixations is 1/mu.  In other words, calling haploid_genome_cleaner every generation
  is an unnecessary run-time cost.  Now, we check if any fixations exist, and
  return
  if there are none.  We check separately for neutral- and non-neutral
  fixations.

  Naively, we could simply scan the mutations container of the entire
  simulation for fixations.
  But, we can do better than that because:

  1. Fixations are, by definition, found in all haploid_genomes.  Thus, it suffices to
  check the
  first (extant) haploid_genome.

  2. When we find the first fixation in a haploid_genome, it is by definition the
  fixation with the smallest
  position, out of all fixations.  We can store its value and use it to look up
  the first fixation
  in all of the remaining haploid_genomes.  Searching for the first fixation in this
  way avoids cache misses that
  are unavoidable when we do out-of-order lookups in "mcounts" for the
  remaining fixations.
*/

namespace fwdpp
{
    namespace fwdpp_internal
    {

        // First, we have a set of helper functions:

        //! Wrapper around std::find_if to find next haploid_genome where n > 0
        template <typename GenomeContainerType_iterator>
        inline GenomeContainerType_iterator
        next_extant_haploid_genome(GenomeContainerType_iterator beg,
                                   GenomeContainerType_iterator end) noexcept
        {
            return std::find_if(
                beg, end,
                [](const typename GenomeContainerType_iterator::value_type &g) noexcept {
                    return g.n;
                });
        }

        struct find_fixation
        //! Find the first mutation in a haploid_genome that is also a fixation
        {
            template <typename mut_index_cont, typename mcounts_t>
            inline typename mut_index_cont::const_iterator
            operator()(const mut_index_cont &mc, const mcounts_t &mcounts,
                       const uint_t twoN) const noexcept
            {
                return std::find_if(mc.cbegin(), mc.cend(),
                                    [&mcounts, &twoN](const std::size_t &i) noexcept {
                                        return mcounts[i] == twoN;
                                    });
            }

            template <typename mut_index_cont, typename MutationContainerType,
                      typename mcounts_t, typename mutation_removal_policy>
            inline typename mut_index_cont::const_iterator
            operator()(const mut_index_cont &mc, const MutationContainerType &mutations,
                       const mcounts_t &mcounts, const uint_t twoN,
                       mutation_removal_policy &&mp) const noexcept
            {
                return std::find_if(
                    mc.cbegin(), mc.cend(),
                    [&mutations, &mcounts, &twoN, &mp](const std::size_t &i) noexcept {
                        return mcounts[i] == twoN && mp(mutations[i]);
                    });
            }
        };

        struct haploid_genome_cleaner_erase_remove_idiom_wrapper
        //! Wrapper for erase/remove idiom.
        {
            template <typename mut_index_cont, typename mcounts_t>
            inline void
            operator()(mut_index_cont &mc, const mcounts_t &mcounts,
                       const typename mut_index_cont::value_type &first_fixation,
                       const uint_t twoN) const noexcept
            //! Overload for no custom policy
            {
                /*
                  The first call to std::find relies on a de-referencing of the
                  return value of
                  find_fixation, passed to this function as first_fixation.
                  Because mc is sorted according to mutation position,
                  first_fixation is the fixation with the smallest position.
                  Using std::find like this avoids some calls to the lambda
                  expression, which experiences cache-misses because mcounts is
                  NOT sorted according to mutation position.
                */
                mc.erase(std::remove_if(
                             std::find(mc.begin(), mc.end(), first_fixation), mc.end(),
                             [&mcounts, &twoN](
                                 const typename mut_index_cont::value_type &i) noexcept {
                                 return mcounts[i] == twoN;
                             }),
                         mc.end());
            }

            template <typename mut_index_cont, typename MutationContainerType,
                      typename mcounts_t, typename mutation_removal_policy>
            inline void
            operator()(mut_index_cont &mc, const MutationContainerType &mutations,
                       const mcounts_t &mcounts,
                       const typename mut_index_cont::value_type &first_fixation,
                       const uint_t twoN, mutation_removal_policy &&mp) const noexcept
            //! Overload for custom policy.
            {
                /*
                  The first call to std::find relies on a de-referencing of the
                  return value of
                  find_fixation, passed to this function as first_fixation.
                  Because mc is sorted according to mutation position,
                  first_fixation is the fixation with the smallest position.
                  Using std::find like this avoids some calls to the lambda
                  expression, which experiences cache-misses because mcounts is
                  NOT sorted according to mutation position.
                */
                mc.erase(std::remove_if(
                             std::find(mc.begin(), mc.end(), first_fixation), mc.end(),
                             [&mcounts, &mutations, &twoN, &mp](
                                 const typename mut_index_cont::value_type &i) noexcept {
                                 return mcounts[i] == twoN && mp(mutations[i]);
                             }),
                         mc.end());
            }

            template <typename mut_index_cont, typename mcounts_t>
            inline void
            operator()(mut_index_cont &mc, const mcounts_t &mcounts,
                       const uint_t twoN) const noexcept
            //! Overload for no custom policy
            {
                /*
                  \note Added in 0.5.0 to address Issue #41
                */
                mc.erase(std::remove_if(
                             mc.begin(), mc.end(),
                             [&mcounts, &twoN](
                                 const typename mut_index_cont::value_type &i) noexcept {
                                 return mcounts[i] == twoN;
                             }),
                         mc.end());
            }

            template <typename mut_index_cont, typename MutationContainerType,
                      typename mcounts_t, typename mutation_removal_policy>
            inline void
            operator()(mut_index_cont &mc, const MutationContainerType &mutations,
                       const mcounts_t &mcounts, const uint_t twoN,
                       mutation_removal_policy &&mp) const noexcept
            //! Overload for custom policy.
            {
                /*
                  \note Added in 0.5.0 to address Issue #41
                */
                mc.erase(std::remove_if(
                             mc.begin(), mc.end(),
                             [&mcounts, &mutations, &twoN, &mp](
                                 const typename mut_index_cont::value_type &i) noexcept {
                                 return mcounts[i] == twoN && mp(mutations[i]);
                             }),
                         mc.end());
            }
        };

        /*! \brief Handles removal of indexes to mutations from haploid_genomes after
          sampling
          Intended use is when std::is_same< mutation_removal_policy,
          fwdpp::remove_nothing >::type is true.
          Called by fwdpp::sample_diploid
        */
        template <typename GenomeContainerType, typename MutationContainerType,
                  typename mutation_removal_policy>
        inline typename std::enable_if<
            std::is_same<mutation_removal_policy, fwdpp::remove_nothing>::value>::type
        haploid_genome_cleaner(GenomeContainerType &, const MutationContainerType &,
                               const std::vector<uint_t> &, const uint_t,
                               const mutation_removal_policy &)
        {
            return;
        }

        template <typename GenomeContainerType, typename fixation_finder,
                  typename idiom_wrapper>
        inline void
        haploid_genome_cleaner_details(GenomeContainerType &haploid_genomes,
                                       const fixation_finder &ff,
                                       const idiom_wrapper &iw) noexcept
        /*!
          The two overloads of haploid_genome_cleaner dispatch the above policies into
          this function.
        */
        {
            auto extant_haploid_genome = next_extant_haploid_genome(
                haploid_genomes.begin(), haploid_genomes.end());
            const auto fixation_n = ff(extant_haploid_genome->mutations);
            bool neutral_fixations_exist
                = (fixation_n != extant_haploid_genome->mutations.cend());
            const auto fixation_s = ff(extant_haploid_genome->smutations);
            bool selected_fixations_exist
                = (fixation_s != extant_haploid_genome->smutations.cend());
            if (!neutral_fixations_exist && !selected_fixations_exist)
                return;

            const auto gend = haploid_genomes.end();
            // Assign values to avoid tons of de-referencing later
            const auto fixation_n_value
                = (fixation_n == extant_haploid_genome->mutations.cend())
                      ? typename decltype(fixation_n)::value_type()
                      : *fixation_n;
            const auto fixation_s_value
                = (fixation_s == extant_haploid_genome->smutations.cend())
                      ? typename decltype(fixation_s)::value_type()
                      : *fixation_s;
            while (extant_haploid_genome < gend)
                {
                    if (neutral_fixations_exist)
                        {
                            iw(extant_haploid_genome->mutations, fixation_n_value);
                        }
                    if (selected_fixations_exist)
                        {
                            iw(extant_haploid_genome->smutations, fixation_s_value);
                        }
                    extant_haploid_genome
                        = next_extant_haploid_genome(extant_haploid_genome + 1, gend);
                }
        }

        template <typename GenomeContainerType, typename fixation_finder>
        std::pair<bool, bool>
        fixation_finder_search_all(const GenomeContainerType &haploid_genomes,
                                   const fixation_finder &ff) noexcept
        /*!
          Determines whether a fixation exists in ANY haploid_genome.
         */
        {
            bool neutral = false, selected = false;
            for (const auto &g : haploid_genomes)
                {
                    if (!neutral)
                        {
                            auto itr = ff(g.mutations);
                            if (itr != g.mutations.cend())
                                {
                                    neutral = true;
                                }
                        }
                    if (!selected)
                        {
                            auto itr = ff(g.smutations);
                            if (itr != g.smutations.cend())
                                {
                                    selected = true;
                                }
                        }
                    if (neutral && selected)
                        return std::make_pair(neutral, selected);
                }
            return std::make_pair(neutral, selected);
        }

        template <typename GenomeContainerType, typename fixation_finder,
                  typename idiom_wrapper>
        inline void
        haploid_genome_cleaner_details(GenomeContainerType &haploid_genomes,
                                       const fixation_finder &ff,
                                       const idiom_wrapper &iw, std::true_type) noexcept
        /*!
          The two overloads of haploid_genome_cleaner dispatch the above policies into
          this function.
          \brief Overload for case where a model must search all haploid_genomes to
          figure out if a fixation exists.
          Use case is the multi-locus API.
          \note Added in 0.5.0 to address issue #41 where the logic of this
          routine failed if the first
          extant haploid_genome was in "locus 1" but the first fixation is in "locus 0"
        */
        {
            auto t = fixation_finder_search_all(haploid_genomes, ff);
            auto neutral_fixations_exist = t.first;
            auto selected_fixations_exist = t.second;
            if (!neutral_fixations_exist && !selected_fixations_exist)
                return;

            auto extant_haploid_genome = next_extant_haploid_genome(
                haploid_genomes.begin(), haploid_genomes.end());
            const auto gend = haploid_genomes.end();
            while (extant_haploid_genome < gend)
                {
                    if (neutral_fixations_exist)
                        {
                            iw(extant_haploid_genome->mutations);
                        }
                    if (selected_fixations_exist)
                        {
                            iw(extant_haploid_genome->smutations);
                        }
                    extant_haploid_genome
                        = next_extant_haploid_genome(extant_haploid_genome + 1, gend);
                }
        }

        /*
          Now, we have the two overloads of haploid_genome_cleaner.
          The std::cref(...) below are required, o/w copies will get bound,
          which is bad for performance.
        */

        /*! \brief Handles removal of indexes to mutations from haploid_genomes after
          sampling
          Intended use is when std::is_same< mutation_removal_policy,
          fwdpp::true_type >::type is true.
          Called by fwdpp::sample_diploid
        */
        template <typename GenomeContainerType, typename MutationContainerType,
                  typename mutation_removal_policy>
        inline typename std::enable_if<
            std::is_same<mutation_removal_policy, std::true_type>::value>::type
        haploid_genome_cleaner(GenomeContainerType &haploid_genomes,
                               const MutationContainerType &,
                               const std::vector<uint_t> &mcounts, const uint_t twoN,
                               const mutation_removal_policy &)
        {
            haploid_genome_cleaner_details(
                haploid_genomes,
                [&mcounts,
                 twoN](const typename GenomeContainerType::value_type::mutation_container
                           &mc) { return find_fixation()(mc, mcounts, twoN); },
                [&mcounts,
                 twoN](typename GenomeContainerType::value_type::mutation_container &mc,
                       const typename GenomeContainerType::value_type::
                           mutation_container::value_type v) {
                    return haploid_genome_cleaner_erase_remove_idiom_wrapper()(
                        mc, mcounts, v, twoN);
                });
        }

        /*! \brief Handles removal of indexes to mutations from haploid_genomes after
          sampling
          This overload handles truly custom policies, which must take a
          mutation type as an argument.
          Called by fwdpp::sample_diploid
        */
        template <typename GenomeContainerType, typename MutationContainerType,
                  typename mutation_removal_policy>
        inline typename std::enable_if<
            !std::is_same<mutation_removal_policy, std::true_type>::value
            && !std::is_same<mutation_removal_policy,
                             fwdpp::remove_nothing>::value>::type
        haploid_genome_cleaner(GenomeContainerType &haploid_genomes,
                               const MutationContainerType &mutations,
                               const std::vector<uint_t> &mcounts, const uint_t twoN,
                               const mutation_removal_policy &mp)
        {
            haploid_genome_cleaner_details(
                haploid_genomes,
                [&mutations, &mcounts, &mp,
                 twoN](const typename GenomeContainerType::value_type::mutation_container
                           &mc) {
                    return find_fixation()(mc, mutations, mcounts, twoN, mp);
                },
                [&mutations, &mcounts, &mp,
                 twoN](typename GenomeContainerType::value_type::mutation_container &mc,
                       typename GenomeContainerType::value_type::mutation_container::
                           value_type v) {
                    return haploid_genome_cleaner_erase_remove_idiom_wrapper()(
                        mc, mutations, mcounts, v, twoN, mp);
                });
        }

        /*! \brief Handles removal of indexes to mutations from haploid_genomes after
          sampling
          Intended use is when std::is_same< mutation_removal_policy,
          fwdpp::true_type >::type is true.
          Called by fwdpp::sample_diploid

          \note Added in 0.5.0 to deal with issue #41 involving multi-locus
          sims
        */
        template <typename GenomeContainerType, typename MutationContainerType,
                  typename mutation_removal_policy>
        inline typename std::enable_if<
            std::is_same<mutation_removal_policy, std::true_type>::value>::type
        haploid_genome_cleaner(GenomeContainerType &haploid_genomes,
                               const MutationContainerType &,
                               const std::vector<uint_t> &mcounts, const uint_t twoN,
                               const mutation_removal_policy &, std::true_type)
        {
            haploid_genome_cleaner_details(
                haploid_genomes,
                [&mcounts,
                 twoN](const typename GenomeContainerType::value_type::mutation_container
                           &mc) { return find_fixation()(mc, mcounts, twoN); },
                [&mcounts, twoN](
                    typename GenomeContainerType::value_type::mutation_container &mc) {
                    return haploid_genome_cleaner_erase_remove_idiom_wrapper()(
                        mc, mcounts, twoN);
                },
                std::true_type());
        }

        /*! \brief Handles removal of indexes to mutations from haploid_genomes after
          sampling
          This overload handles truly custom policies, which must take a
          mutation type as an argument.
          Called by fwdpp::sample_diploid

          \note Added in 0.5.0 to deal with issue #41 involving multi-locus
          sims
        */
        template <typename GenomeContainerType, typename MutationContainerType,
                  typename mutation_removal_policy>
        inline typename std::enable_if<
            !std::is_same<mutation_removal_policy, std::true_type>::value
            && !std::is_same<mutation_removal_policy,
                             fwdpp::remove_nothing>::value>::type
        haploid_genome_cleaner(GenomeContainerType &haploid_genomes,
                               const MutationContainerType &mutations,
                               const std::vector<uint_t> &mcounts, const uint_t twoN,
                               const mutation_removal_policy &mp, std::true_type)
        {
            haploid_genome_cleaner_details(
                haploid_genomes,
                [&mutations, &mcounts, &mp,
                 twoN](const typename GenomeContainerType::value_type::mutation_container
                           &mc) {
                    return find_fixation()(mc, mutations, mcounts, twoN, mp);
                },
                [&mutations, &mcounts, &mp, twoN](
                    typename GenomeContainerType::value_type::mutation_container &mc) {
                    return haploid_genome_cleaner_erase_remove_idiom_wrapper()(
                        mc, mutations, mcounts, twoN, mp);
                },
                std::true_type());
        }
    }
}

#endif
