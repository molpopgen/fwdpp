#ifndef FWDPP_ALGORITHM_COMPACT_MUTATIONS_HPP__
#define FWDPP_ALGORITHM_COMPACT_MUTATIONS_HPP__

#include <vector>
#include <algorithm>
#include <numeric>
#include <fwdpp/forward_types.hpp>

namespace fwdpp
{
    template <typename T>
    std::vector<fwdpp::uint_t>
    compact_mutations(T &pop)
    /// \brief Sort mutation positions in population.
    ///
    /// Reorders the population mutation container
    /// so that it is sorted by increasing mutation position.
    /// The reordering requires assigning new key values into
    /// all haploid_genomes, which is also done. The mutation counts
    /// and mut_lookup data structures also get updated.
    ///
    /// Running this periodically is a performance increase
    /// due to improved cache performance.
    ///
    /// \param pop A population.
    ///
    /// \return A vector containing the mapping of old to
    /// new mutation index
    {
        std::vector<fwdpp::uint_t> indexes(pop.mutations.size());
        std::iota(std::begin(indexes), std::end(indexes), 0);

        auto new_indexes_end = std::stable_partition(
            std::begin(indexes), std::end(indexes),
            [&pop](const fwdpp::uint_t i) { return pop.mcounts[i]; });

        std::sort(std::begin(indexes), new_indexes_end,
                  [&pop](const fwdpp::uint_t i, const fwdpp::uint_t j) {
                      return pop.mutations[i].pos < pop.mutations[j].pos;
                  });
        std::vector<fwdpp::uint_t> reindex(indexes.size());
        std::size_t new_indexes_size
            = std::distance(std::begin(indexes), new_indexes_end);
        for (std::size_t i = 0; i < new_indexes_size; ++i)
            {
                reindex[indexes[i]] = i;
            }
        for (auto &g : pop.haploid_genomes)
            {
                if (g.n)
                    {
                        for (auto &m : g.mutations)
                            {
                                m = reindex[m];
                            }
                        for (auto &m : g.smutations)
                            {
                                m = reindex[m];
                            }
                    }
            }
        decltype(pop.mutations) reordered_muts;
        decltype(pop.mcounts) reordered_mcounts;
        reordered_muts.reserve(pop.mutations.size());
        reordered_mcounts.reserve(pop.mutations.size());
        for (auto i : indexes)
            {
                reordered_muts.emplace_back(std::move(pop.mutations[i]));
                reordered_mcounts.push_back(pop.mcounts[i]);
                if (reordered_mcounts.back() > 0)
                    {
                        auto x = pop.mut_lookup.equal_range(
                            reordered_muts.back().pos);
                        while (x.first != x.second)
                            {
                                if (x.first->second == i)
                                    {
                                        x.first->second
                                            = reordered_muts.size() - 1;
                                        break;
                                    }
                                ++x.first;
                            }
                    }
            }
        pop.mutations.swap(reordered_muts);
        pop.mcounts.swap(reordered_mcounts);
        return reindex;
    }
} // namespace fwdpp

#endif
