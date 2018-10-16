#ifndef FWDPP_TS_MARGINAL_TREE_ITERATOR_HPP
#define FWDPP_TS_MARGINAL_TREE_ITERATOR_HPP

#include "marginal_tree.hpp"
#include "table_collection.hpp"
#include "iterate_marginal_trees_details.hpp"

namespace fwdpp
{
    namespace ts
    {
        class marginal_tree_iterator
        {
          private:
            indexed_edge_container::const_iterator j, jM, k, kM;
            double x, maxpos;

          public:
            marginal_tree marginal;
            marginal_tree_iterator(const table_collection& tables,
                                   const std::vector<TS_NODE_INT>& samples)
                : j(tables.input_left.cbegin()), jM(tables.input_left.cend()),
                  k(tables.output_right.cbegin()),
                  kM(tables.output_right.cend()), x(0.0), maxpos(tables.L),
                  marginal(tables.num_nodes(), samples)
            {
            }

            marginal_tree_iterator(
                const table_collection& tables,
                const std::vector<TS_NODE_INT>& samples,
                const std::vector<TS_NODE_INT>& preserved_nodes)
                : j(tables.input_left.cbegin()), jM(tables.input_left.cend()),
                  k(tables.output_right.cbegin()),
                  kM(tables.output_right.cend()), x(0.0), maxpos(tables.L),
                  marginal(tables.num_nodes(), samples, preserved_nodes)
            {
            }

            template <typename leaf_policy, typename sample_list_policy>
            inline bool
            operator()(const leaf_policy lp, const sample_list_policy slp)
            {
                if (j < jM || x < maxpos)
                    {
                        while (k < kM && k->pos == x) // T4
                            {
                                const auto p = k->parent;
                                const auto c = k->child;
                                const auto lsib = marginal.left_sib[c];
                                const auto rsib = marginal.right_sib[c];
                                if (lsib == TS_NULL_NODE)
                                    {
                                        marginal.left_child[p] = rsib;
                                    }
                                else
                                    {
                                        marginal.right_sib[lsib] = rsib;
                                    }
                                if (rsib == TS_NULL_NODE)
                                    {
                                        marginal.right_child[p] = lsib;
                                    }
                                else
                                    {
                                        marginal.left_sib[rsib] = lsib;
                                    }
                                marginal.parents[c] = TS_NULL_NODE;
                                marginal.left_sib[c] = TS_NULL_NODE;
                                marginal.right_sib[c] = TS_NULL_NODE;
                                detail::outgoing_leaf_counts(
                                    marginal, k->parent, k->child, lp);
                                detail::update_sample_list(marginal, k->parent,
                                                           slp);
                                ++k;
                            }
                        while (j < jM && j->pos == x) // Step T2
                            {
                                const auto p = j->parent;
                                const auto c = j->child;
                                const auto rchild = marginal.right_child[p];
                                if (rchild == TS_NULL_NODE)
                                    {
                                        marginal.left_child[p] = c;
                                        marginal.left_sib[c] = TS_NULL_NODE;
                                        marginal.right_sib[c] = TS_NULL_NODE;
                                    }
                                else
                                    {
                                        marginal.right_sib[rchild] = c;
                                        marginal.left_sib[c] = rchild;
                                        marginal.right_sib[c] = TS_NULL_NODE;
                                    }
                                // The entry for the child refers to
                                // the parent's location in the node table.
                                marginal.parents[c] = j->parent;
                                marginal.right_child[p] = c;
                                detail::incoming_leaf_counts(
                                    marginal, j->parent, j->child, lp);
                                detail::update_sample_list(marginal, j->parent,
                                                           slp);
                                ++j;
                            }

                        double right = maxpos;
                        if (j < jM)
                            {
                                right = std::min(right, j->pos);
                            }
                        if (k < kM)
                            {
                                right = std::min(right, k->pos);
                            }
                        marginal.left = x;
                        marginal.right = right;
                        bool rv = j < jM || x < maxpos;
                        x = right;
                        return rv;
                    }
                return false;
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
