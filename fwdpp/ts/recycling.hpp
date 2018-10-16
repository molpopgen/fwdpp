#ifndef FWDPP_TS_RECYCLING_HPP
#define FWDPP_TS_RECYCLING_HPP

#include <vector>
#include <cstdint>
#include <queue>

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
    } // namespace ts
} // namespace fwdpp

#endif
