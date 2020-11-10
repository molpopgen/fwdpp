/// \file util.hpp
#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/gsl_discrete.hpp>
#include <fwdpp/internal/util.hpp>
#include <set>
#include <map>
#include <type_traits>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace fwdpp
{
    template <typename poptype>
    void
    zero_out_haploid_genomes(poptype &pop)
    /// Set haploid_genome counts in all diploids to zero
    /// \version 0.7.4 Added to fwdpp
    {
        fwdpp_internal::zero_out_haploid_genomes(pop, typename poptype::popmodel_t());
    }

    /*!
      Label all extinct and fixed variants for recycling

      \note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename MutationContainerType, typename mutation_lookup_table>
    void
    update_mutations(MutationContainerType &mutations, mutation_lookup_table &lookup,
                     std::vector<uint_t> &mcounts, const unsigned twoN)
    {
        static_assert(traits::is_mutation_v<typename MutationContainerType::value_type>,
                      "mutation_type must be derived from fwdpp::mutation_base");
#ifndef NDEBUG
        if (mcounts.size() != mutations.size())
            {
                throw std::runtime_error("FWDPP DEBUG: mutation counts size "
                                         "must equal mutation container size");
            }
#endif
        for (std::size_t i = 0; i < mcounts.size(); ++i)
            {
#ifndef NDEBUG
                if (mcounts[i] > twoN)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: mutation count out of range");
                    }
#endif
                if (mcounts[i] == twoN || !mcounts[i])
                    {
                        auto itr = lookup.equal_range(mutations[i].pos);
                        while (itr.first != itr.second)
                            {
                                if (itr.first->second == i)
                                    {
                                        lookup.erase(itr.first);
                                        mcounts[i] = 0;
                                        break;
                                    }
                                ++itr.first;
                            }
                    }
            }
    }

    /*!
      Label all extinct  variants for recycling

      \note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename MutationContainerType, typename mutation_lookup_table>
    void
    update_mutations(const MutationContainerType &mutations,
                     mutation_lookup_table &lookup, std::vector<uint_t> &mcounts)

    {
        static_assert(traits::is_mutation_v<typename MutationContainerType::value_type>,
                      "mutation_type must be derived from fwdpp::mutation_base");
        for (std::size_t i = 0; i < mcounts.size(); ++i)
            {
                if (!mcounts[i])
                    {
                        auto itr = lookup.equal_range(mutations[i].pos);
                        while (itr.first != itr.second)
                            {
                                if (itr.first->second == i)
                                    {
                                        lookup.erase(mutations[i].pos);
                                        mcounts[i] = 0;
                                        break;
                                    }
                                ++itr.first;
                            }
                    }
            }
    }

    /*!
      Label all fixed and all extinct variants for recycling. Copy fixations
      and fixation times
      into containers.

      \note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename MutationContainerType, typename fixation_container_t,
              typename fixation_time_container_t, typename mutation_lookup_table>
    void
    update_mutations(MutationContainerType &mutations, fixation_container_t &fixations,
                     fixation_time_container_t &fixation_times,
                     mutation_lookup_table &lookup, std::vector<uint_t> &mcounts,
                     const unsigned &generation, const unsigned &twoN)
    {
        static_assert(traits::is_mutation_v<typename MutationContainerType::value_type>,
                      "mutation_type must be derived from fwdpp::mutation_base");
#ifndef NDEBUG
        if (mcounts.size() != mutations.size())
            {
                throw std::runtime_error("FWDPP DEBUG: mutation counts size "
                                         "must equal mutation container size");
            }
#endif
        for (unsigned i = 0; i < mcounts.size(); ++i)
            {
#ifndef NDEBUG
                if (mcounts[i] > twoN)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: mutation count out of range");
                    }
#endif
                if (mcounts[i] == twoN)
                    {
                        fixations.push_back(mutations[i]);
                        fixation_times.push_back(generation);
                        mcounts[i] = 0; // set count to zero to mark mutation
                                        // as "recyclable"
                        auto itr = lookup.equal_range(mutations[i].pos);
                        while (itr.first != itr.second)
                            {
                                if (itr.first->second == i)
                                    {
                                        lookup.erase(itr.first);
                                        break;
                                    }
                                ++itr.first;
                            }
                    }
                else if (!mcounts[i])
                    {
                        auto itr = lookup.equal_range(mutations[i].pos);
                        if (itr.first != lookup.end())
                            {
                                while (itr.first != itr.second)
                                    {
                                        if (itr.first->second == i)
                                            {
                                                lookup.erase(itr.first);
                                                break;
                                            }
                                        ++itr.first;
                                    }
                            }
                    }
            }
    }

    /*!
      Label all fixed neutral variant and all extinct variants for recycling.
      Copy fixations and fixation times
      for neutral mutations into containers.

      \note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename MutationContainerType, typename fixation_container_t,
              typename fixation_time_container_t, typename mutation_lookup_table>
    void
    update_mutations_n(MutationContainerType &mutations, fixation_container_t &fixations,
                       fixation_time_container_t &fixation_times,
                       mutation_lookup_table &lookup, std::vector<uint_t> &mcounts,
                       const unsigned &generation, const unsigned &twoN)
    {
        static_assert(traits::is_mutation_v<typename MutationContainerType::value_type>,
                      "mutation_type must be derived from "
                      "fwdpp::mutation_base");
#ifndef NDEBUG
        if (mcounts.size() != mutations.size())
            {
                throw std::runtime_error("FWDPP DEBUG: mutation counts size "
                                         "must equal mutation container size");
            }
#endif
        for (unsigned i = 0; i < mcounts.size(); ++i)
            {
#ifndef NDEBUG
                if (mcounts[i] > twoN)
                    {
                        throw std::runtime_error(
                            "FWDPP DEBUG: mutation count out of range");
                    }
#endif
                if (mutations[i].neutral && mcounts[i] == twoN)
                    {
                        fixations.push_back(mutations[i]);
                        fixation_times.push_back(generation);
                        mcounts[i] = 0; // set count to zero to mark mutation
                                        // as "recyclable"
                        auto itr = lookup.equal_range(mutations[i].pos);
                        while (itr.first != itr.second)
                            {
                                if (itr.first->second == i)
                                    {
                                        lookup.erase(mutations[i].pos);
                                        break;
                                    }
                                ++itr.first;
                            }
                    }
                else if (mcounts[i] == twoN) // Fixation is not neutral
                    {
                        // Guard against repeatedly recording fixation data
                        if (std::find(std::begin(fixations), std::end(fixations),
                                      mutations[i])
                            == std::end(fixations))
                            {
                                fixations.push_back(mutations[i]);
                                fixation_times.push_back(generation);
                            }
                    }
                else if (!mcounts[i])
                    {
                        auto itr = lookup.equal_range(mutations[i].pos);
                        while (itr.first != itr.second)
                            {
                                if (itr.first->second == i)
                                    {
                                        lookup.erase(mutations[i].pos);
                                        break;
                                    }
                                ++itr.first;
                            }
                    }
            }
    }
} // namespace fwdpp
#endif /* _UTIL_HPP_ */
