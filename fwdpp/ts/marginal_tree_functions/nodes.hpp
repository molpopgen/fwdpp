#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODES_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODES_HPP

#include <vector>
#include <algorithm>
#include <memory>
#include "roots.hpp"
#include "children.hpp"
#include "node_traversal_order.hpp"
#include "node_traversal_preorder.hpp"
#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        class node_iterator
        /// \brief Traverse nodes in a marginal_tree
        /// \headerfile fwdpp/ts/marginal_tree_functions/nodes.hpp
        {
          private:
            const marginal_tree &t;
            std::vector<table_index_t> subtree_roots;
            table_index_t current_root;
            std::unique_ptr<node_traversal_order> order;

            std::vector<table_index_t>
            init_subtree_roots()
            {
                auto r = get_roots(t);
                std::reverse(begin(r), end(r));
                return r;
            }

            std::vector<table_index_t>
            init_subtree_roots(table_index_t u)
            {
                if (static_cast<std::size_t>(u) >= t.left_child.size())
                    {
                        throw std::invalid_argument("node it out of range");
                    }
                return { u };
            }

            table_index_t
            init_current_root()
            {
                if (subtree_roots.empty())
                    {
                        return NULL_INDEX;
                    }
                auto rv = subtree_roots.back();
                subtree_roots.pop_back();
                return rv;
            }

          public:
            template <typename ORDER>
            node_iterator(const marginal_tree &m, ORDER order_policy)
                : t(m), subtree_roots(init_subtree_roots()),
                  current_root(init_current_root()),
                  order(node_traversal_dispatch(current_root, order_policy))
            /// Traverse nodes starting from all roots
            /// \param m A marginal_tree
            /// \param order_policy A dispatch tag.
            ///
            /// The dispatch tag is sent to the appropriate overload of
            /// fwdpp::ts::node_traversal_dispatch, which must return a std::unique_ptr
            /// containing an object derived from node_traversal_order.
            {
            }

            template <typename ORDER>
            node_iterator(const marginal_tree &m, table_index_t u,
                          ORDER order_policy)
                : t(m), subtree_roots(init_subtree_roots(u)),
                  current_root(init_current_root()),
                  order(node_traversal_dispatch(current_root, order_policy))
            /// Traverse nodes in a subtree whose root is \a u
            /// \param m A marginal_tree
            /// \param u A node in \a m
            /// \param order_policy A dispatch tag.
            ///
            /// The dispatch tag is sent to the appropriate overload of
            /// fwdpp::ts::node_traversal_dispatch, which must return a std::unique_ptr
            /// containing an object derived from node_traversal_order.
            {
            }

            inline table_index_t
            operator()()
            /// Return the next node in the tree.
            /// A value ot NULL_INDEX signals end of iteration.
            {
                if (current_root != NULL_INDEX)
                    {
                        auto rv = order->operator()(t);
                        if (rv == NULL_INDEX)
                            {
                                current_root = init_current_root();
                                if (current_root != NULL_INDEX)
                                    {
                                        order->initialize(current_root);
                                        rv = order->operator()(t);
                                    }
                                else
                                    {
                                        rv = current_root;
                                    }
                            }
                        return rv;
                    }
                return NULL_INDEX;
            }

            template <typename F>
            inline table_index_t
            operator()(const F &f)
            /// Apply a function to each node.
            /// \param f A function behaving as void(*process_node)(table_index_t)
            /// Returns false to signal end of iteration.
            {
                auto v = this->operator()();
                bool rv = (v != NULL_INDEX);
                if (rv)
                    {
                        f(v);
                    }
                return rv;
            }
        };

        template <typename ORDER, typename F>
        inline void
        process_nodes(const marginal_tree &m, ORDER order, const F &f)
        /// Apply a function to all nodes in the tree \a m.
        /// \param m A marginal_tree
        /// \param order A dispatch tag specifying the node traversal roder
        /// \param f A function behaving as void(*process_node)(table_index_t)
        {
            node_iterator ni(m, order);
            while (ni(f))
                {
                }
        }

        template <typename ORDER, typename F>
        inline void
        process_nodes(const marginal_tree &m, table_index_t u, ORDER order,
                      const F &f)
        /// Apply a function to all nodes in the subtree whose root is \a u
        /// \param m A marginal_tree
        /// \param u A node id in \a m
        /// \param order A dispatch tag specifying the node traversal roder
        /// \param f A function behaving as void(*process_node)(table_index_t)
        {
            node_iterator ni(m, u, order);
            while (ni(f))
                {
                }
        }

        inline int
        num_nodes(const marginal_tree &m)
        /// Get the number of nodes in \a m
        /// \param m A marginal_tree
        {
            int nnodes = 0;
            process_nodes(m, nodes_preorder(),
                          [&nnodes](const table_index_t /*n*/) { ++nnodes; });
            return nnodes;
        }

        template <typename ORDER>
        inline std::vector<table_index_t>
        get_nodes(const marginal_tree &m, ORDER order)
        /// Get a vector of all nodes in the tree \a m
        /// \param m A marginal_tree
        /// \param order A dispatch tag specifying the node traversal roder
        {
            std::vector<table_index_t> nodes;
            process_nodes(m, order, [&nodes](const table_index_t n) {
                nodes.push_back(n);
            });
            return nodes;
        }

        template <typename ORDER>
        inline std::vector<table_index_t>
        get_nodes(const marginal_tree &m, table_index_t u, ORDER order)
        /// Get a vector of all nodes in the subtree of \a m whose root is \a u
        /// \param m A marginal_tree
        /// \param u A node in \a m
        /// \param order A dispatch tag specifying the node traversal roder
        {
            std::vector<table_index_t> nodes;
            process_nodes(m, u, order, [&nodes](const table_index_t n) {
                nodes.push_back(n);
            });
            return nodes;
        }

    } // namespace ts
} // namespace fwdpp

#endif
