/// \file util.hpp
#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <set>
#include <map>
#include <type_traits>
#include <algorithm>
#include <functional>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace KTfwd
{
    /*!
      Label all extinct and fixed variants for recycling

      \note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename mcont_t, typename mutation_lookup_table>
    void
    update_mutations(mcont_t &mutations, mutation_lookup_table &lookup,
                     std::vector<uint_t> &mcounts, const unsigned twoN)
    {
        static_assert(
            typename traits::is_mutation<typename mcont_t::value_type>::type(),
            "mutation_type must be derived from KTfwd::mutation_base");
        assert(mcounts.size() == mutations.size());
        for (std::size_t i = 0; i < mcounts.size(); ++i)
            {
                assert(mcounts[i] <= twoN);
                if (mcounts[i] == twoN || !mcounts[i])
                    {
                        lookup.erase(mutations[i].pos);
                        mcounts[i] = 0;
                    }
            }
    }

    /*!
      Label all extinct  variants for recycling

      \Note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename mcont_t, typename mutation_lookup_table>
    void
    update_mutations(const mcont_t &mutations, mutation_lookup_table &lookup,
                     std::vector<uint_t> &mcounts)

    {
        static_assert(
            typename traits::is_mutation<typename mcont_t::value_type>::type(),
            "mutation_type must be derived from KTfwd::mutation_base");
        for (std::size_t i = 0; i < mcounts.size(); ++i)
            {
                if (!mcounts[i])
                    {
                        lookup.erase(mutations[i].pos);
                        mcounts[i] = 0;
                    }
            }
    }

    /*!
      Label all fixed and all extinct variants for recycling. Copy fixations
      and fixation times
      into containers.

      \note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename mcont_t, typename fixation_container_t,
              typename fixation_time_container_t,
              typename mutation_lookup_table>
    void
    update_mutations(mcont_t &mutations, fixation_container_t &fixations,
                     fixation_time_container_t &fixation_times,
                     mutation_lookup_table &lookup,
                     std::vector<uint_t> &mcounts, const unsigned &generation,
                     const unsigned &twoN)
    {
        static_assert(
            typename traits::is_mutation<typename mcont_t::value_type>::type(),
            "mutation_type must be derived from KTfwd::mutation_base");
        assert(mcounts.size() == mutations.size());
        for (unsigned i = 0; i < mcounts.size(); ++i)
            {
                assert(mcounts[i] <= twoN);
                if (mcounts[i] == twoN)
                    {
                        fixations.push_back(mutations[i]);
                        fixation_times.push_back(generation);
                        mcounts[i] = 0; // set count to zero to mark mutation
                                        // as "recyclable"
                        lookup.erase(mutations[i].pos);
                    }
                if (!mcounts[i])
                    lookup.erase(mutations[i].pos);
            }
    }

    /*!
      Label all fixed neutral variant and all extinct variants for recycling.
      Copy fixations and fixation times
      for neutral mutations into containers.

      \note: lookup must be compatible with lookup->erase(lookup->find(double))
    */
    template <typename mcont_t, typename fixation_container_t,
              typename fixation_time_container_t,
              typename mutation_lookup_table>
    void
    update_mutations_n(mcont_t &mutations, fixation_container_t &fixations,
                       fixation_time_container_t &fixation_times,
                       mutation_lookup_table &lookup,
                       std::vector<uint_t> &mcounts,
                       const unsigned &generation, const unsigned &twoN)
    {
        static_assert(
            typename traits::is_mutation<typename mcont_t::value_type>::type(),
            "mutation_type must be derived from KTfwd::mutation_base");
        assert(mcounts.size() == mutations.size());
        for (unsigned i = 0; i < mcounts.size(); ++i)
            {
                assert(mcounts[i] <= twoN);
                if (mutations[i].neutral && mcounts[i] == twoN)
                    {
                        fixations.push_back(mutations[i]);
                        fixation_times.push_back(generation);
                        mcounts[i] = 0; // set count to zero to mark mutation
                                        // as "recyclable"
                        lookup.erase(mutations[i].pos);
                    }
                if (!mcounts[i])
                    lookup.erase(mutations[i].pos);
            }
    }
}
#endif /* _UTIL_HPP_ */
