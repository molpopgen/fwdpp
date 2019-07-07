#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODE_TRAVERSAL_ORDER_HPP__
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_NODE_TRAVERSAL_ORDER_HPP__

#include "../marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        struct node_traversal_order
        /// \brief Interface class for dependency injection
        ///        into node_iterator
        /// \headerfile fwdpp/ts/marginal_tree_functions/node_traversal_order.hpp
        {
            node_traversal_order() = default;
            virtual ~node_traversal_order() = default;
            node_traversal_order(const node_traversal_order&) = default;
            node_traversal_order(node_traversal_order&&) = default;
            node_traversal_order& operator=(const node_traversal_order&)
                = default;
            node_traversal_order& operator=(node_traversal_order&&) = default;

            virtual TS_NODE_INT operator()(const marginal_tree&) = 0;
            virtual void initialize(TS_NODE_INT /*root*/) = 0;
        };
    } // namespace ts
} // namespace fwdpp

#endif
