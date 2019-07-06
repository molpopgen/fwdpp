#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODES_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODES_HPP

#include <stack>
#include <vector>
#include "roots.hpp"
#include "children.hpp"
#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            template <typename STACK>
            inline TS_NODE_INT
            node_traversal_preorder(const marginal_tree &m, STACK &nstack)
            {
                if (!nstack.empty())
                    {
                        auto rv = nstack.top();
                        nstack.pop();
                        if (num_children(m, rv) != 0)
                            {
                                // traverse children right to left
                                process_children(m, rv, false,
                                                 [&nstack](TS_NODE_INT x) {
                                                     nstack.push(x);
                                                 });
                            }
                        return rv;
                    }
                return TS_NULL_NODE;
            }
        }; // namespace detail

        struct nodes_preorder
        {
        };

        struct nodes_inorder
        {
        };

        template <typename STACK>
        inline TS_NODE_INT
        node_traversal_dispatch(const marginal_tree &m, STACK &nstack,
                                nodes_preorder)
        {
            return detail::node_traversal_preorder(m, nstack);
        }

        class node_iterator
        {
          private:
            // NOTE: default is std::deque, but there is no way that that
            // is the optimal choice.
            using node_stack
                = std::stack<TS_NODE_INT, std::vector<TS_NODE_INT>>;

            const marginal_tree &t;
            root_iterator ri;
            TS_NODE_INT current_root, current_node;
            node_stack nstack;

          public:
            explicit node_iterator(const marginal_tree &m)
                : t(m), ri(m), current_root(ri()),
                  current_node(current_root), nstack{ { current_root } }
            {
            }

            template <typename ORDER>
            inline TS_NODE_INT
            operator()(ORDER order)
            {
                if (current_root != TS_NULL_NODE)
                    {
                        current_node
                            = node_traversal_dispatch(t, nstack, order);
                        if (current_node == TS_NULL_NODE)
                            {
                                current_root = ri();
                                return current_root;
                            }
                        return current_node;
                    }
                current_root = ri();
                return current_root; // must be TS_NULL_NODE
            }

            template <typename ORDER, typename F>
            inline TS_NODE_INT
            operator()(ORDER order, const F &f)
            {
                auto v = this->operator()(order);
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
            node_iterator ni(m);
            while (ni(order, f))
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

    } // namespace ts
} // namespace fwdpp

#endif
