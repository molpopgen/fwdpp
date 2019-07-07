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
        {
          private:
            const marginal_tree &t;
            std::vector<TS_NODE_INT> subtree_roots;
            TS_NODE_INT current_root;
            std::unique_ptr<node_traversal_order> order;

            std::vector<TS_NODE_INT>
            init_subtree_roots()
            {
                auto r = get_roots(t);
                std::reverse(begin(r), end(r));
                return r;
            }

            std::vector<TS_NODE_INT>
            init_subtree_roots(TS_NODE_INT u)
            {
                if (static_cast<std::size_t>(u) >= t.left_child.size())
                    {
                        throw std::invalid_argument("node it out of range");
                    }
                return { u };
            }

            TS_NODE_INT
            init_current_root()
            {
                if (subtree_roots.empty())
                    {
                        return TS_NULL_NODE;
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
            {
            }

            template <typename ORDER>
            node_iterator(const marginal_tree &m, TS_NODE_INT u,
                          ORDER order_policy)
                : t(m), subtree_roots(init_subtree_roots(u)),
                  current_root(init_current_root()),
                  order(node_traversal_dispatch(current_root, order_policy))
            {
            }

            inline TS_NODE_INT
            operator()()
            {
                if (current_root != TS_NULL_NODE)
                    {
                        auto rv = order->operator()(t);
                        if (rv == TS_NULL_NODE)
                            {
                                current_root = init_current_root();
                                if (current_root != TS_NULL_NODE)
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
                return TS_NULL_NODE;
            }

            template <typename F>
            inline TS_NODE_INT
            operator()(const F &f)
            {
                auto v = this->operator()();
                bool rv = (v != TS_NULL_NODE);
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
        {
            node_iterator ni(m, order);
            while (ni(f))
                {
                }
        }

        template <typename ORDER, typename F>
        inline void
        process_nodes(const marginal_tree &m, TS_NODE_INT u, ORDER order,
                      const F &f)
        {
            node_iterator ni(m, u, order);
            while (ni(f))
                {
                }
        }

        inline int
        num_nodes(const marginal_tree &m)
        {
            int nnodes = 0;
            process_nodes(m, nodes_preorder(),
                          [&nnodes](const TS_NODE_INT /*n*/) { ++nnodes; });
            return nnodes;
        }

        template <typename ORDER>
        inline std::vector<TS_NODE_INT>
        get_nodes(const marginal_tree &m, ORDER order)
        {
            std::vector<TS_NODE_INT> nodes;
            process_nodes(m, order, [&nodes](const TS_NODE_INT n) {
                nodes.push_back(n);
            });
            return nodes;
        }

        template <typename ORDER>
        inline std::vector<TS_NODE_INT>
        get_nodes(const marginal_tree &m, TS_NODE_INT u, ORDER order)
        {
            std::vector<TS_NODE_INT> nodes;
            process_nodes(m, u, order, [&nodes](const TS_NODE_INT n) {
                nodes.push_back(n);
            });
            return nodes;
        }

    } // namespace ts
} // namespace fwdpp

#endif
