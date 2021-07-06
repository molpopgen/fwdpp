#ifndef FWDPP_TS_TREE_VISITOR_HPP
#define FWDPP_TS_TREE_VISITOR_HPP

#include <vector>
#include <algorithm>
#include "exceptions.hpp"
#include "marginal_tree.hpp"
#include "detail/advance_marginal_tree_policies.hpp"
#include <fwdpp/util/named_type.hpp>

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

        template <typename TableCollectionType> class tree_visitor
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
            using id_type = typename TableCollectionType::id_type;
            static constexpr id_type null = TableCollectionType::null;
            typename std::vector<id_type>::const_iterator j, jM, k, kM;
            typename TableCollectionType::edge_table::const_iterator beg_edges,
                end_edges;
            double x, maxpos;
            marginal_tree<id_type> marginal;
            bool advancing_sample_list;

            void
            update_roots_outgoing(id_type p, id_type c, marginal_tree<id_type>& marginal)
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
                        while (x != null && above_sample == 0)
                            {
                                // If this node is a sample, then
                                // it must descend from a root,
                                // which is why we OR this in the loop below.
                                above_sample = (marginal.sample_index_map[x] != null);
                                auto lc = marginal.left_child[x];
                                while (lc != null && above_sample == 0)
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
                                marginal.left_root = null;
                                if (lroot != null)
                                    {
                                        marginal.right_sib[lroot] = rroot;
                                        marginal.left_root = lroot;
                                    }
                                if (rroot != null)
                                    {
                                        marginal.left_sib[rroot] = lroot;
                                        marginal.left_root = rroot;
                                    }
                                marginal.left_sib[root] = null;
                                marginal.right_sib[root] = null;
                            }
                        if (marginal.left_root != null)
                            {
                                //Put c into a pre-existing root list
                                auto lroot = marginal.left_sib[marginal.left_root];
                                if (lroot != null)
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
            update_roots_incoming(id_type p, id_type c, id_type lsib, id_type rsib,
                                  marginal_tree<id_type>& marginal)
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
                        while (x != null && above_sample == 0)
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
                                if (lsib != null)
                                    {
                                        marginal.right_sib[lsib] = root;
                                    }
                                if (rsib != null)
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
                                marginal.left_root = null;
                                if (lsib != null)
                                    {
                                        marginal.right_sib[lsib] = rsib;
                                        marginal.left_root = lsib;
                                    }
                                if (rsib != null)
                                    {
                                        marginal.left_sib[rsib] = lsib;
                                        marginal.left_root = rsib;
                                    }
                            }
                    }
            }

          public:
            template <typename SAMPLES>
            tree_visitor(const TableCollectionType& tables, SAMPLES&& samples,
                         update_samples_list update)
                : j(tables.input_left.cbegin()), jM(tables.input_left.cend()),
                  k(tables.output_right.cbegin()), kM(tables.output_right.cend()),
                  beg_edges(begin(tables.edges)), end_edges(end(tables.edges)), x(0.0),
                  maxpos(tables.genome_length()),
                  marginal(tables.num_nodes(), std::forward<SAMPLES>(samples),
                           update.get()),
                  advancing_sample_list(update.get())
            /// \todo Document
            {
                if ((j == jM || k == kM) && !tables.edges.empty())
                    {
                        throw std::invalid_argument("tables are not indexed");
                    }
            }

            tree_visitor(const TableCollectionType& tables,
                         const std::vector<id_type>& samples,
                         const std::vector<id_type>& preserved_nodes,
                         update_samples_list update)
                : j(tables.input_left.cbegin()), jM(tables.input_left.cend()),
                  k(tables.output_right.cbegin()), kM(tables.output_right.cend()),
                  beg_edges(begin(tables.edges)), end_edges(end(tables.edges)), x(0.0),
                  maxpos(tables.genome_length()),
                  marginal(tables.num_nodes(), samples, preserved_nodes, update.get()),
                  advancing_sample_list(update.get())
            {
                if ((j == jM || k == kM) && !tables.edges.empty())
                    {
                        throw std::invalid_argument("tables are not indexed");
                    }
                if (samples.empty() && preserved_nodes.empty())
                    {
                        throw samples_error("one or both sample lists are empty");
                    }
            }

            const marginal_tree<typename TableCollectionType::id_type>&
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

            inline bool
            operator()()
            {
                if (j < jM || x < maxpos)
                    {
                        while (k < kM && (beg_edges + *k)->right == x) // T4
                            {
                                const auto p = (beg_edges + *k)->parent;
                                const auto c = (beg_edges + *k)->child;
                                const auto lsib = marginal.left_sib[c];
                                const auto rsib = marginal.right_sib[c];
                                if (lsib == null)
                                    {
                                        marginal.left_child[p] = rsib;
                                    }
                                else
                                    {
                                        marginal.right_sib[lsib] = rsib;
                                    }
                                if (rsib == null)
                                    {
                                        marginal.right_child[p] = lsib;
                                    }
                                else
                                    {
                                        marginal.left_sib[rsib] = lsib;
                                    }
                                marginal.parents[c] = null;
                                marginal.left_sib[c] = null;
                                marginal.right_sib[c] = null;
                                detail::outgoing_leaf_counts(marginal,
                                                             (beg_edges + *k)->parent,
                                                             (beg_edges + *k)->child);
                                if (advancing_sample_list)
                                    {
                                        detail::update_samples_list(
                                            marginal, (beg_edges + *k)->parent);
                                    }
                                update_roots_outgoing(p, c, marginal);
                                ++k;
                            }
                        while (j < jM && (beg_edges + *j)->left == x) // Step T2
                            {
                                const auto p = (beg_edges + *j)->parent;
                                const auto c = (beg_edges + *j)->child;
                                const auto rchild = marginal.right_child[p];
                                const auto lsib = marginal.left_sib[c];
                                const auto rsib = marginal.right_sib[c];
                                if (rchild == null)
                                    {
                                        marginal.left_child[p] = c;
                                        marginal.left_sib[c] = null;
                                        marginal.right_sib[c] = null;
                                    }
                                else
                                    {
                                        marginal.right_sib[rchild] = c;
                                        marginal.left_sib[c] = rchild;
                                        marginal.right_sib[c] = null;
                                    }
                                // The entry for the child refers to
                                // the parent's location in the node table.
                                marginal.parents[c] = (beg_edges + *j)->parent;
                                marginal.right_child[p] = c;
                                detail::incoming_leaf_counts(marginal,
                                                             (beg_edges + *j)->parent,
                                                             (beg_edges + *j)->child);
                                if (advancing_sample_list)
                                    {
                                        detail::update_samples_list(
                                            marginal, (beg_edges + *j)->parent);
                                    }
                                update_roots_incoming(p, c, lsib, rsib, marginal);

                                ++j;
                            }

                        // This is a big "gotcha".
                        // The root tracking functions will sometimes
                        // result in left_root not actually being the left_root.
                        // We loop through the left_sibs to fix that.
                        if (marginal.left_root != null)
                            {
                                while (marginal.left_sib[marginal.left_root] != null)
                                    {
                                        marginal.left_root
                                            = marginal.left_sib[marginal.left_root];
                                    }
                            }
#ifndef NDEBUG
                        // Validate the roots via brute-force.
                        auto lr = marginal.left_root;
                        if (lr == null)
                            {
                                throw std::runtime_error(
                                    "FWDPP DEBUG: left_root is null");
                            }
                        std::vector<int> is_root(marginal.sample_index_map.size(), 0);
                        std::vector<int> processed(is_root.size(), 0);
                        for (std::size_t s = 0; s < marginal.sample_index_map.size();
                             ++s)
                            {
                                if (marginal.sample_index_map[s] != null)
                                    {
                                        auto u = static_cast<id_type>(s);
                                        auto root = u;
                                        bool early_exit = false;
                                        while (u != null)
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
                                throw std::runtime_error("FWDPP DEBUG: num_roots "
                                                         "disagreement");
                            }
                        while (lr != null)
                            {
                                if (is_root[lr] != 1)
                                    {
                                        throw std::runtime_error("FWDPP DEBUG: root "
                                                                 "contents "
                                                                 "disagreement");
                                    }
                                lr = marginal.right_sib[lr];
                            }
#endif
                        double right = maxpos;
                        if (j < jM)
                            {
                                right = std::min(right, (beg_edges + *j)->left);
                            }
                        if (k < kM)
                            {
                                right = std::min(right, (beg_edges + *k)->right);
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
