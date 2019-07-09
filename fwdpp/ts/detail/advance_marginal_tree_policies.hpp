#ifndef FWDPP_TS_ADVANCE_MARGINAL_TREE_POLICIES_HPP
#define FWDPP_TS_ADVANCE_MARGINAL_TREE_POLICIES_HPP

#include <cstdint>
#include <type_traits>
#include <cassert>
#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            inline void
            outgoing_leaf_counts(marginal_tree& marginal,
                                 const std::int32_t parent,
                                 const std::int32_t child)
            {
                auto p = parent;
                auto lc = marginal.leaf_counts[child];
                auto plc = marginal.preserved_leaf_counts[child];
                while (p != TS_NULL_NODE)
                    {
                        marginal.leaf_counts[p] -= lc;
                        marginal.preserved_leaf_counts[p] -= plc;
                        assert(marginal.leaf_counts[p] >= 0);
                        assert(marginal.preserved_leaf_counts[p] >= 0);
                        p = marginal.parents[p];
                    }
            }

            inline void
            incoming_leaf_counts(marginal_tree& marginal,
                                 const std::int32_t parent,
                                 const std::int32_t child)
            {
                auto p = parent;
                auto lc = marginal.leaf_counts[child];
                auto plc = marginal.preserved_leaf_counts[child];
                while (p != TS_NULL_NODE)
                    {
                        marginal.leaf_counts[p] += lc;
                        marginal.preserved_leaf_counts[p] += plc;
                        p = marginal.parents[p];
                    }
            }

            inline void
            update_samples_list(marginal_tree& marginal,
                                const std::int32_t node)
            {
                const auto& parents = marginal.parents;
                const auto& sample_map = marginal.sample_index_map;
                const auto& left_child = marginal.left_child;
                const auto& right_sib = marginal.right_sib;

                auto& right = marginal.right_sample;
                auto& left = marginal.left_sample;
                auto& next = marginal.next_sample;
                for (auto n = node; n != TS_NULL_NODE; n = parents[n])
                    {
                        auto sample_index = sample_map[n];
                        if (sample_index != TS_NULL_NODE)
                            {
                                right[n] = left[n];
                            }
                        else
                            {
                                left[n] = TS_NULL_NODE;
                                right[n] = TS_NULL_NODE;
                            }
                        for (auto v = left_child[n]; v != TS_NULL_NODE;
                             v = right_sib[v])
                            {
                                if (left[v] != TS_NULL_NODE)
                                    {
                                        assert(right[v] != TS_NULL_NODE);
                                        if (left[n] == TS_NULL_NODE)
                                            {
                                                left[n] = left[v];
                                                right[n] = right[v];
                                            }
                                        else
                                            {
                                                next[right[n]] = left[v];
                                                right[n] = right[v];
                                            }
                                    }
                            }
                    }
            }
        } // namespace detail
    }     // namespace ts
} // namespace fwdpp

#endif
