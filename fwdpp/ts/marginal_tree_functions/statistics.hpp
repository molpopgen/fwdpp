#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_STATISTICS_HPP
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_STATISTICS_HPP

#include <stdexcept>
#include "../marginal_tree.hpp"
#include "../node.hpp"
#include "nodes.hpp"

namespace fwdpp
{
    namespace ts
    {
        double
        total_time(const marginal_tree &m, const std::vector<node> &nodes,
                   bool scale_by_length)
        {
            if (m.parents.size() != nodes.size())
                {
                    throw std::invalid_argument(
                        "inconsistent container sizes");
                }
            double ttime = 0.;
            process_nodes(
                m, nodes_preorder(), [&ttime, &m, &nodes](table_index_t u) {
                    if (m.parents[u] != NULL_INDEX)
                        {
                            ttime
                                += (nodes[u].time - nodes[m.parents[u]].time);
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
