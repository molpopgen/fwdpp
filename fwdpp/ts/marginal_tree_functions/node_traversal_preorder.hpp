#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODE_TRAVERSAL_PREORDER_HPP__
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODE_TRAVERSAL_PREORDER_HPP__

#include <stack>
#include <vector>
#include <memory>
#include "children.hpp"
#include "node_traversal_order.hpp"

namespace fwdpp
{
    namespace ts
    {
        class node_traversal_preorder : public node_traversal_order
        /// \brief Preorder traversal of nodes for a node_iterator
        /// \headerfile fwdpp/ts/marginal_tree_functions/node_traversal_preorder.hpp
        {
          private:
            using node_stack
                = std::stack<table_index_t, std::vector<table_index_t>>;
            node_stack nstack;
            table_index_t current_node;

          public:
            explicit node_traversal_preorder(table_index_t u)
                : node_traversal_order(), nstack{ { u } }, current_node{
                      NULL_INDEX
                  }
            {
            }

            table_index_t
            operator()(const marginal_tree& m) final
            {
                if (!nstack.empty())
                    {
                        current_node = nstack.top();
                        nstack.pop();
                        if (num_children(m, current_node) != 0)
                            {
                                process_children(
                                    m, current_node, false,
                                    [this](table_index_t x) { nstack.push(x); });
                            }
                        return current_node;
                    }
                return NULL_INDEX;
            }

            void
            initialize(table_index_t root) final
            {
                nstack.push(root);
            }
        };

        struct nodes_preorder
        /// Dispatch tage for node_traversal_preorder
        {
        };

        inline std::unique_ptr<node_traversal_order>
        node_traversal_dispatch(table_index_t root, nodes_preorder)
        /// Handles dependency injection of node_traversal_preorder into node_iterator
        {
            return std::unique_ptr<node_traversal_order>(
                new node_traversal_preorder(root));
        }

    } // namespace ts
} // namespace fwdpp

#endif

