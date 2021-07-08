#ifndef FWDPP_TS_ADVANCE_MARGINAL_TREE_POLICIES_HPP
#define FWDPP_TS_ADVANCE_MARGINAL_TREE_POLICIES_HPP

#include <cstdint>
#include <type_traits>
#include <cassert>
#include "../types/generate_null_id.hpp"
#include "../marginal_tree.hpp"
#include "../types/generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            template <typename SignedInteger>
            void
            update_roots_outgoing(SignedInteger p, SignedInteger c,
                                  marginal_tree<SignedInteger>& marginal)
            // This is the algorithm used by tskit.
            // The method is applied to nodes p and c
            // AFTER the output index has been applied
            // to a marginal tree.
            //
            // The criteria to mark something as a root
            // are that its parent is NULL and it is
            // above a sample. Here, c's parent is NULL
            // because that is the result of processing
            // "outgoing" nodes.  Thus, c will be
            // added to the list of nodes.
            //
            // The trick is to determine if the most
            // ancient ancestor of p is NOT a root,
            // which occurs if that node is not above any samples.
            //
            // While we are at it, we update above_sample.
            // A node is above a sample if it is a sample
            // or if any of its children are above a sample.
            {
                if (marginal.above_sample[c] == 1)
                    {
                        auto x = p;
                        auto root = x;
                        std::int8_t above_sample = 0;
                        while (x != types::generate_null_id<SignedInteger>()
                               && above_sample == 0)
                            {
                                // If this node is a sample, then
                                // it must descend from a root,
                                // which is why we OR this in the loop below.
                                above_sample
                                    = (marginal.sample_index_map[x]
                                       != types::generate_null_id<SignedInteger>());
                                auto lc = marginal.left_child[x];
                                while (lc != types::generate_null_id<SignedInteger>()
                                       && above_sample == 0)
                                    {
                                        above_sample
                                            = above_sample || marginal.above_sample[lc];
                                        lc = marginal.right_sib[lc];
                                    }
                                marginal.above_sample[x] = above_sample;
                                root = x;
                                x = marginal.parents[x];
                            }
                        // Now, root refers to the most ancient
                        // ancestor of p found in the above loop
                        if (above_sample == 0)
                            {
                                // Remove root from list of roots
                                auto lroot = marginal.left_sib[root];
                                auto rroot = marginal.right_sib[root];
                                marginal.left_root
                                    = types::generate_null_id<SignedInteger>();
                                if (lroot != types::generate_null_id<SignedInteger>())
                                    {
                                        marginal.right_sib[lroot] = rroot;
                                        marginal.left_root = lroot;
                                    }
                                if (rroot != types::generate_null_id<SignedInteger>())
                                    {
                                        marginal.left_sib[rroot] = lroot;
                                        marginal.left_root = rroot;
                                    }
                                marginal.left_sib[root]
                                    = types::generate_null_id<SignedInteger>();
                                marginal.right_sib[root]
                                    = types::generate_null_id<SignedInteger>();
                            }
                        if (marginal.left_root
                            != types::generate_null_id<SignedInteger>())
                            {
                                //Put c into a pre-existing root list
                                auto lroot = marginal.left_sib[marginal.left_root];
                                if (lroot != types::generate_null_id<SignedInteger>())
                                    {
                                        marginal.right_sib[lroot] = c;
                                    }
                                marginal.left_sib[c] = lroot;
                                marginal.left_sib[marginal.left_root] = c;
                            }
                        marginal.right_sib[c] = marginal.left_root;
                        marginal.left_root = c;
                    }
            }

            template <typename SignedInteger>
            void
            update_roots_incoming(SignedInteger p, SignedInteger c, SignedInteger lsib,
                                  SignedInteger rsib,
                                  marginal_tree<SignedInteger>& marginal)
            // This is the algorithm used by tskit.  It is applied
            // AFTER the incoming node list has updated marginal.
            // Thus, p is parent to c.
            //
            // lsib and rsib are with respect to c BEFORE the incoming
            // node list has updated marginal.  This is a bit confusing:
            // these values are used to remove c from the root list if
            // necessary.
            //
            // NOTE: all values of c for which above_sample[c] == 1
            // are in the root list.
            {
                if (marginal.above_sample[c])
                    {
                        auto x = p;
                        auto root = x;

                        std::int8_t above_sample = 0;
                        while (x != types::generate_null_id<SignedInteger>()
                               && above_sample == 0)
                            {
                                above_sample = marginal.above_sample[x];
                                // c is above_sample and p is c's parent.
                                // Thus, all parents to p are above_sample, too.
                                marginal.above_sample[x] = 1;
                                root = x;
                                x = marginal.parents[x];
                            }
                        if (above_sample == 0)
                            // If we are here, then the above loop terminated
                            // by encountering a null node, because above_sample[x]
                            // must have been 0 for all x. However, because c is
                            // above sample, all nodes encountered have been update
                            // to be above_sample as well. Thus, the new value of root
                            // replaces c in the root list.
                            {
                                // Replace c with root in root list
                                if (lsib != types::generate_null_id<SignedInteger>())
                                    {
                                        marginal.right_sib[lsib] = root;
                                    }
                                if (rsib != types::generate_null_id<SignedInteger>())
                                    {
                                        marginal.left_sib[rsib] = root;
                                    }
                                marginal.left_sib[root] = lsib;
                                marginal.right_sib[root] = rsib;
                                marginal.left_root = root;
                            }
                        else
                            // If we are here, then we encountered a node
                            // ancestral to c where above_sample == 1.
                            // Thus, c can no longer be a root.  If the current
                            // p is also a c in a later call to this function, then
                            // it may also be removed, etc..
                            {
                                // Remove c from root list
                                marginal.left_root
                                    = types::generate_null_id<SignedInteger>();
                                if (lsib != types::generate_null_id<SignedInteger>())
                                    {
                                        marginal.right_sib[lsib] = rsib;
                                        marginal.left_root = lsib;
                                    }
                                if (rsib != types::generate_null_id<SignedInteger>())
                                    {
                                        marginal.left_sib[rsib] = lsib;
                                        marginal.left_root = rsib;
                                    }
                            }
                    }
            }

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
