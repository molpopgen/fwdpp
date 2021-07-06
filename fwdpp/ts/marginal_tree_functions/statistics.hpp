#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_STATISTICS_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_STATISTICS_HPP

#include <stdexcept>
#include "../marginal_tree.hpp"
#include "../types/generate_null_id.hpp"
#include "../types/node.hpp"
#include "nodes.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger>
        double
        total_time(const marginal_tree<SignedInteger> &m,
                   const std::vector<types::node<SignedInteger>> &nodes,
                   bool scale_by_length)
        {
            if (m.parents.size() != nodes.size())
                {
                    throw std::invalid_argument("inconsistent container sizes");
                }
            double ttime = 0.;
            process_nodes(m, nodes_preorder(), [&ttime, &m, &nodes](SignedInteger u) {
                if (m.parents[u] != types::generate_null_id<SignedInteger>())
                    {
                        ttime += (nodes[u].time - nodes[m.parents[u]].time);
                    }
            });
            if (scale_by_length)
                {
                    ttime *= (m.right - m.left);
                }
            return ttime;
        }
    } // namespace ts
} // namespace fwdpp

#endif
