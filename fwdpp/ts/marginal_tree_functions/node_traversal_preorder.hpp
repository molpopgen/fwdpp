#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODE_TRAVERSAL_PREORDER_HPP__
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODE_TRAVERSAL_PREORDER_HPP__

#include <stack>
#include <vector>
#include <memory>
#include "children.hpp"
#include "node_traversal_order.hpp"
#include "../types/generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger>
        class node_traversal_preorder : public node_traversal_order<SignedInteger>
        /// \brief Preorder traversal of nodes for a node_iterator
        /// \headerfile fwdpp/ts/marginal_tree_functions/node_traversal_preorder.hpp
        {
          private:
            using node_stack = std::stack<SignedInteger, std::vector<SignedInteger>>;
            node_stack nstack;
            SignedInteger current_node;

          public:
            explicit node_traversal_preorder(SignedInteger u)
                : node_traversal_order<SignedInteger>(), nstack{{u}},
                  current_node{types::generate_null_id<SignedInteger>()}
            {
            }

            SignedInteger
            operator()(const marginal_tree<SignedInteger>& m) final
            {
                if (!nstack.empty())
                    {
                        current_node = nstack.top();
                        nstack.pop();
                        if (num_children(m, current_node) != 0)
                            {
                                process_children(
                                    m, current_node, false,
                                    [this](SignedInteger x) { nstack.push(x); });
                            }
                        return current_node;
                    }
                return types::generate_null_id<SignedInteger>();
            }

            void
            initialize(SignedInteger root) final
            {
                nstack.push(root);
            }
        };

        struct nodes_preorder
        /// Dispatch tage for node_traversal_preorder
        {
        };

        template <typename SignedInteger>
        inline std::unique_ptr<node_traversal_order<SignedInteger>>
        node_traversal_dispatch(SignedInteger root, nodes_preorder)
        /// Handles dependency injection of node_traversal_preorder into node_iterator
        {
            return std::unique_ptr<node_traversal_order<SignedInteger>>(
                new node_traversal_preorder<SignedInteger>(root));
        }

    } // namespace ts
} // namespace fwdpp

#endif

