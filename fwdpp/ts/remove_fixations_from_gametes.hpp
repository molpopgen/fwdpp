#ifndef FWDPP_TS_REMOVE_FIXATIONS_FROM_GAMETES_HPP
#define FWDPP_TS_REMOVE_FIXATIONS_FROM_GAMETES_HPP

#include <algorithm>
#include <cstdint>
#include <fwdpp/forward_types.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename gcont_t, typename mcont_t,
                  typename mutation_count_container>
        void
        remove_fixations_from_gametes(
            gcont_t &gametes, const mcont_t &mutations,
            const mutation_count_container &mcounts,
            const mutation_count_container &mcounts_from_preserved_nodes,
            const fwdpp::uint_t twoN)
        {
            bool fixations_exist = false;
            for (std::size_t i = 0; !fixations_exist && i < mcounts.size();
                 ++i)
                {
                    if (mcounts[i] == twoN
                        && mcounts_from_preserved_nodes[i] == 0)
                        {
                            fixations_exist = true;
                        }
                }
            if (fixations_exist)
                {
                    auto removal_criteria
                        = [&mcounts, &mcounts_from_preserved_nodes,
                           twoN](const fwdpp::uint_t key) {
                              return mcounts[key] == twoN
                                     && mcounts_from_preserved_nodes[key] == 0;
                          };
                    for (auto &g : gametes)
                        {
                            if (g.n)
                                {
                                    auto itr = std::remove_if(
                                        g.mutations.begin(), g.mutations.end(),
                                        removal_criteria);
                                    g.mutations.erase(itr, g.mutations.end());
                                    itr = std::remove_if(g.smutations.begin(),
                                                         g.smutations.end(),
                                                         removal_criteria);
                                    g.smutations.erase(itr,
                                                       g.smutations.end());
                                }
                        }
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
