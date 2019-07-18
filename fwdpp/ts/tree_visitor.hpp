#ifndef FWDPP_TS_TREE_VISITOR_HPP
#define FWDPP_TS_TREE_VISITOR_HPP

#include <vector>
#include <algorithm>
#include "exceptions.hpp"
#include "marginal_tree.hpp"
#include "table_collection.hpp"
#include "detail/advance_marginal_tree_policies.hpp"
#include <fwdpp/named_type.hpp>

namespace fwdpp
{
    namespace ts
    {
        struct update_samples_list_t
        {
        };

        /// Policy dictating if a tree_visitor updates sample lists.
        using update_samples_list
            = strong_types::named_type<bool, update_samples_list_t>;

        class tree_visitor
        /// \brief Class that iterates over marginal trees.
        ///
        /// \note This class declares private data
        /// whose integrity are tied to the lifetime
        /// of the table_collection used to construct
        /// a tree_visitor!
        ///
        /// \version 0.7.0 Added to fwdpp
        /// \version 0.7.4 Updates tree roots during traversal.
        {
          private:
            indexed_edge_container::const_iterator j, jM, k, kM;
            double x, maxpos;
            marginal_tree marginal;
            bool advancing_sample_list;

            void
            update_roots_outgoing(TS_NODE_INT p, TS_NODE_INT c,
                                  marginal_tree& marginal)
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
                        while (x != TS_NULL_NODE && above_sample == 0)
                            {
                                // If this node is a sample, then
                                // it must descend from a root,
                                // which is why we OR this in the loop below.
                                above_sample = (marginal.sample_index_map[x]
                                                != TS_NULL_NODE);
                                auto lc = marginal.left_child[x];
                                while (lc != TS_NULL_NODE && above_sample == 0)
                                    {
                                        above_sample
                                            = above_sample
                                              || marginal.above_sample[lc];
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
                                marginal.left_root = TS_NULL_NODE;
                                if (lroot != TS_NULL_NODE)
                                    {
                                        marginal.right_sib[lroot] = rroot;
                                        marginal.left_root = lroot;
                                    }
                                if (rroot != TS_NULL_NODE)
                                    {
                                        marginal.left_sib[rroot] = lroot;
                                        marginal.left_root = rroot;
                                    }
                                marginal.left_sib[root] = TS_NULL_NODE;
                                marginal.right_sib[root] = TS_NULL_NODE;
                            }
                        if (marginal.left_root != TS_NULL_NODE)
                            {
                                //Put c into a pre-existing root list
                                auto lroot
                                    = marginal.left_sib[marginal.left_root];
                                if (lroot != TS_NULL_NODE)
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

            void
            update_roots_incoming(TS_NODE_INT p, TS_NODE_INT c,
                                  TS_NODE_INT lsib, TS_NODE_INT rsib,
                                  marginal_tree& marginal)
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
                        while (x != TS_NULL_NODE && above_sample == 0)
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
                            // by encountering a NULL node, because above_sample[x]
                            // must have been 0 for all x. However, because c is
                            // above sample, all nodes encountered have been update
                            // to be above_sample as well. Thus, the new value of root
                            // replaces c in the root list.
                            {
                                // Replace c with root in root list
                                if (lsib != TS_NULL_NODE)
                                    {
                                        marginal.right_sib[lsib] = root;
                                    }
                                if (rsib != TS_NULL_NODE)
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
                                marginal.left_root = TS_NULL_NODE;
                                if (lsib != TS_NULL_NODE)
                                    {
                                        marginal.right_sib[lsib] = rsib;
                                        marginal.left_root = lsib;
                                    }
                                if (rsib != TS_NULL_NODE)
                                    {
                                        marginal.left_sib[rsib] = lsib;
                                        marginal.left_root = rsib;
                                    }
                            }
                    }
            }

          public:
            template <typename SAMPLES>
            tree_visitor(const table_collection& tables, SAMPLES&& samples,
                         update_samples_list update)
                : j(tables.input_left.cbegin()), jM(tables.input_left.cend()),
                  k(tables.output_right.cbegin()),
                  kM(tables.output_right.cend()), x(0.0),
                  maxpos(tables.genome_length()),
                  marginal(tables.num_nodes(), std::forward<SAMPLES>(samples),
                           update.get()),
                  advancing_sample_list(update.get())
            /// \todo Document
            {
                if ((j == jM || k == kM) && !tables.edge_table.empty())
                    {
                        throw std::invalid_argument("tables are not indexed");
                    }
            }

            const marginal_tree&
            tree() const
            /*! \brief Returns a handle to the current tree.
             *
			 * \return const reference to the current tree.
			 *
			 * \code{cpp}
			 * // Copy-free "view" of
			 * // the stored fwdpp::ts::marginal_tree
			 * auto & tree = tv.tree();
			 * \endcode
			 */
            {
                return marginal;
            }

            tree_visitor(const table_collection& tables,
                         const std::vector<TS_NODE_INT>& samples,
                         const std::vector<TS_NODE_INT>& preserved_nodes,
                         update_samples_list update)
                : j(tables.input_left.cbegin()), jM(tables.input_left.cend()),
                  k(tables.output_right.cbegin()),
                  kM(tables.output_right.cend()), x(0.0),
                  maxpos(tables.genome_length()),
                  marginal(tables.num_nodes(), samples, preserved_nodes,
                           update.get()),
                  advancing_sample_list(update.get())
            {
                if ((j == jM || k == kM) && !tables.edge_table.empty())
                    {
                        throw std::invalid_argument("tables are not indexed");
                    }
                if (samples.empty() && preserved_nodes.empty())
                    {
                        throw samples_error(
                            "one or both sample lists are empty");
                    }
            }

            inline bool
            operator()()
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
                                    marginal, k->parent, k->child);
                                if (advancing_sample_list)
                                    {
                                        detail::update_samples_list(marginal,
                                                                    k->parent);
                                    }
                                update_roots_outgoing(p, c, marginal);
                                ++k;
                            }
                        while (j < jM && j->pos == x) // Step T2
                            {
                                const auto p = j->parent;
                                const auto c = j->child;
                                const auto rchild = marginal.right_child[p];
                                const auto lsib = marginal.left_sib[c];
                                const auto rsib = marginal.right_sib[c];
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
                                    marginal, j->parent, j->child);
                                if (advancing_sample_list)
                                    {
                                        detail::update_samples_list(marginal,
                                                                    j->parent);
                                    }
                                update_roots_incoming(p, c, lsib, rsib,
                                                      marginal);

                                ++j;
                            }

                        // This is a big "gotcha".
                        // The root tracking functions will sometimes
                        // result in left_root not actually being the left_root.
                        // We loop through the left_sibs to fix that.
                        if (marginal.left_root != TS_NULL_NODE)
                            {
                                while (marginal.left_sib[marginal.left_root]
                                       != TS_NULL_NODE)
                                    {
                                        marginal.left_root
                                            = marginal.left_sib
                                                  [marginal.left_root];
                                    }
                            }
#ifndef NDEBUG
                        // Validate the roots via brute-force.
                        auto lr = marginal.left_root;
                        if (lr == TS_NULL_NODE)
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: left_root is null");
                            }
                        std::vector<int> is_root(
                            marginal.sample_index_map.size(), 0);
                        std::vector<int> processed(is_root.size(), 0);
                        for (std::size_t s = 0;
                             s < marginal.sample_index_map.size(); ++s)
                            {
                                if (marginal.sample_index_map[s]
                                    != TS_NULL_NODE)
                                    {
                                        TS_NODE_INT u = s;
                                        auto root = u;
                                        bool early_exit = false;
                                        while (u != TS_NULL_NODE)
                                            {
                                                if (processed[u])
                                                    {
                                                        early_exit = true;
                                                        break;
                                                    }
                                                processed[u] = 1;
                                                root = u;
                                                u = marginal.parents[u];
                                            }
                                        if (early_exit == false)
                                            {
                                                is_root[root] = 1;
                                            }
                                    }
                            }
                        int nroots_brute = 0;
                        for (auto r : is_root)
                            {
                                nroots_brute += r;
                            }
                        if (nroots_brute != marginal.num_roots())
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: num_roots "
                                    "disagreement");
                            }
                        while (lr != TS_NULL_NODE)
                            {
                                if (is_root[lr] != 1)
                                    {
                                        throw std::runtime_error(
                                            "FWDPP DEBUG: root "
                                            "contents "
                                            "disagreement");
                                    }
                                lr = marginal.right_sib[lr];
                            }
#endif
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
                        // Must set return value before
                        // updating right, else the
                        // last tree will be skipped.
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
