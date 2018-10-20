#ifndef FWDPP_TS_RECYCLING_HPP
#define FWDPP_TS_RECYCLING_HPP

#include <stdexcept>
#include <vector>
#include <cstdint>
#include <queue>
#include <fwdpp/forward_types.hpp>

namespace fwdpp
{
    namespace ts
    {
        std::queue<std::size_t>
        make_mut_queue(
            const std::vector<std::uint32_t> &mcounts,
            const std::vector<std::uint32_t> &counts_from_preserved_nodes)
        /// \brief Make a mutation recycling queue for simulations with tree sequences
        /// \param mcounts Contribution of extant nodes to mutation counts
        /// \param counts_from_preserved_nodes Contribution of extinct nodes to mutation counts
        ///
        /// \returns std::queue<std::size_t>
        {
            std::queue<std::size_t> mutation_recycling_bin;
            for (std::size_t i = 0; i < mcounts.size(); ++i)
                {
                    if (mcounts[i] + counts_from_preserved_nodes[i] == 0)
                        {
                            mutation_recycling_bin.push(i);
                        }
                }
            return mutation_recycling_bin;
        }

        template <typename mcont_t, typename lookup_table,
                  typename mutation_count_container>
        void
        flag_mutations_for_recycling(
            const mcont_t &mutations, mutation_count_container &mcounts,
            mutation_count_container &mcounts_from_preserved_nodes,
            lookup_table &lookup, const fwdpp::uint_t twoN)
        {
            for (std::size_t i = 0; i < mcounts.size(); ++i)
                {
                    if (mcounts_from_preserved_nodes[i] == 0)
                        {
                            if (mcounts[i] > twoN)
                                {
                                    throw std::runtime_error(
                                        "mutation count out of range");
                                }
                            if (mcounts[i] == twoN)
                                {
                                    auto itr
                                        = lookup.equal_range(mutations[i].pos);
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
                            else if (mcounts[i] == 0)
                                {
                                    auto itr
                                        = lookup.equal_range(mutations[i].pos);
                                    if (itr.first != lookup.end())
                                        {
                                            while (itr.first != itr.second)
                                                {
                                                    if (itr.first->second == i)
                                                        {
                                                            lookup.erase(
                                                                itr.first);
                                                            break;
                                                        }
                                                    ++itr.first;
                                                }
                                        }
                                }
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
