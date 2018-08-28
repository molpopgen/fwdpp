#ifndef __SAMPLING_FUNCTIONS_HPP__
#define __SAMPLING_FUNCTIONS_HPP__

#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <utility>
#include <string>
#include <cassert>
#include <type_traits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/data_matrix.hpp>

/*! @defgroup samplingPops Functions related to taking samples from simulated
  populations
*/

namespace fwdpp
{
    inline void
    remove_fixed_keys(std::vector<std::pair<std::size_t, uint_t>> &keys,
                      const uint_t nsam)
    /// Removes all keys were key.second == nsam.
    {
        const auto f = [nsam](const std::pair<std::size_t, uint_t> &p) {
            return p.second == nsam;
        };
        keys.erase(std::remove_if(keys.begin(), keys.end(), f), keys.end());
    }

    template <typename mcont_t>
    inline void
    sort_keys(const mcont_t &mutations,
              std::vector<std::pair<std::size_t, uint_t>> &keys)
    /// Sorts keys by position
    {
        const auto comp
            = [&mutations](const std::pair<std::size_t, uint_t> &a,
                           const std::pair<std::size_t, uint_t> &b) {
                  return mutations[a.first].pos < mutations[b.first].pos;
              };
        std::sort(keys.begin(), keys.end(), comp);
    }

    template <typename poptype>
    auto
    generate_filter_sort_keys(const poptype &pop,
                              const std::vector<std::size_t> &individuals,
                              const bool include_neutral,
                              const bool include_selected,
                              const bool remove_fixed)
        -> decltype(mutation_keys(pop, individuals, include_neutral,
                                  include_selected))
    {
        auto keys = mutation_keys(pop, individuals, include_neutral,
                                  include_selected);
        if (remove_fixed)
            {
                remove_fixed_keys(keys.first, 2 * individuals.size());
                remove_fixed_keys(keys.second, 2 * individuals.size());
            }
        sort_keys(pop.mutations, keys.first);
        sort_keys(pop.mutations, keys.second);
        return keys;
    }

    template <typename poptype>
    data_matrix
    sample_individuals(const poptype &pop,
                       const std::vector<std::size_t> &individuals,
                       const bool include_neutral, const bool include_selected,
                       const bool remove_fixed)
    {
        auto keys = generate_filter_sort_keys(
            pop, individuals, include_neutral, include_selected, remove_fixed);
        return haplotype_matrix(pop, individuals, keys.first, keys.second);
    }

    template <typename poptype>
    std::vector<data_matrix>
    sample_individuals_by_window(const poptype &pop,
                                const std::vector<std::size_t> &individuals,
                                const std::vector<std::pair<double,double>> window_boundaries,
                                const bool include_neutral,
                                const bool include_selected,
                                const bool remove_fixed)
    {
        auto keys = generate_filter_sort_keys(
            pop, individuals, include_neutral, include_selected, remove_fixed);
        std::vector<data_matrix> rv;
        decltype(keys.first.begin()) nstart, nend, sstart, send;

        const auto lbf = [&pop](const std::pair<std::size_t, uint_t> &p,
                                const double pos) {
            return pop.mutations[p.first].pos < pos;
        };

        const auto ubf = [&pop](const double pos,
                                const std::pair<std::size_t, uint_t> &p) {
            return pos < pop.mutations[p.first].pos;
        };
        decltype(keys.first) nwk, swk;
        for (const auto &b : window_boundaries)
            {
                nstart = std::lower_bound(keys.first.begin(), keys.first.end(),
                                          b.first, lbf);
                nend = std::upper_bound(nstart, keys.first.end(), b.second,
                                        ubf);
                sstart = std::lower_bound(keys.second.begin(),
                                          keys.first.end(), b.first, lbf);
                send = std::upper_bound(sstart, keys.second.end(), b.second,
                                        ubf);
                nwk.assign(nstart, nend);
                swk.assign(sstart, send);

                rv.emplace_back(haplotype_matrix(pop, individuals, keys.first,
                                                 keys.second));
            }
        return rv;
    }
}
#endif
