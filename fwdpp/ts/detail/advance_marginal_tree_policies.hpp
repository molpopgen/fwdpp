#ifndef FWDPP_TS_ADVANCE_MARGINAL_TREE_POLICIES_HPP
#define FWDPP_TS_ADVANCE_MARGINAL_TREE_POLICIES_HPP

#include <cstdint>
#include <type_traits>
#include <cassert>
#include "../types/generate_null_id.hpp"
#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            template <typename SignedInteger>
            inline void
            outgoing_leaf_counts(marginal_tree<SignedInteger>& marginal,
                                 const SignedInteger parent, const SignedInteger child)
            {
                auto p = parent;
                auto lc = marginal.leaf_counts[child];
                auto plc = marginal.preserved_leaf_counts[child];
                if (lc + plc == 0)
                    {
                        return;
                    }
                while (p != types::generate_null_id<SignedInteger>())
                    {
                        marginal.leaf_counts[p] -= lc;
                        marginal.preserved_leaf_counts[p] -= plc;
                        assert(marginal.leaf_counts[p] >= 0);
                        assert(marginal.preserved_leaf_counts[p] >= 0);
                        p = marginal.parents[p];
                    }
            }

            template <typename SignedInteger>
            inline void
            incoming_leaf_counts(marginal_tree<SignedInteger>& marginal,
                                 const SignedInteger parent, const SignedInteger child)
            {
                auto p = parent;
                auto lc = marginal.leaf_counts[child];
                auto plc = marginal.preserved_leaf_counts[child];
                if (lc + plc == 0)
                    {
                        return;
                    }
                while (p != types::generate_null_id<SignedInteger>())
                    {
                        marginal.leaf_counts[p] += lc;
                        marginal.preserved_leaf_counts[p] += plc;
                        p = marginal.parents[p];
                    }
            }

            template <typename SignedInteger>
            inline void
            update_samples_list(marginal_tree<SignedInteger>& marginal,
                                const SignedInteger node)
            {
                const auto& parents = marginal.parents;
                const auto& sample_map = marginal.sample_index_map;
                const auto& left_child = marginal.left_child;
                const auto& right_sib = marginal.right_sib;

                auto& right = marginal.right_sample;
                auto& left = marginal.left_sample;
                auto& next = marginal.next_sample;
                for (auto n = node; n != types::generate_null_id<SignedInteger>();
                     n = parents[n])
                    {
                        auto sample_index = sample_map[n];
                        if (sample_index != types::generate_null_id<SignedInteger>())
                            {
                                right[n] = left[n];
                            }
                        else
                            {
                                left[n] = types::generate_null_id<SignedInteger>();
                                right[n] = types::generate_null_id<SignedInteger>();
                            }
                        for (auto v = left_child[n];
                             v != types::generate_null_id<SignedInteger>();
                             v = right_sib[v])
                            {
                                if (left[v] != types::generate_null_id<SignedInteger>())
                                    {
                                        assert(
                                            right[v]
                                            != types::generate_null_id<SignedInteger>());
                                        if (left[n]
                                            == types::generate_null_id<SignedInteger>())
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
