#ifndef FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP
#define FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP

#include <map>
#include <utility>
#include <vector>
#include "definitions.hpp"
#include "table_collection.hpp"
#include "tree_visitor.hpp"

namespace fwdpp
{
    namespace ts
    {
        // TODO: consider flattening the return value to a vector
        inline std::map<TS_NODE_INT, std::vector<std::pair<double, double>>>
        mark_multiple_roots(const table_collection &tables,
                            const std::vector<TS_NODE_INT> &samples)
        /// \brief Identify root nodes in "marginal forests".
        ///
        /// \version 0.7.0 Added to library
        /// \version 0.7.4 Refactored to use root tracking method
        ///
        /// See fwdpp::ts::mutate_tables for discussion.
        {
            std::map<TS_NODE_INT, std::vector<std::pair<double, double>>> rv;
            tree_visitor mti(tables, samples);
            while (mti(std::true_type(), std::false_type()))
                {
                    auto &tree = mti.tree();
                    if (tree.num_roots() > 1)
                        {
                            auto lr = tree.left_root;
                            while (lr != TS_NULL_NODE)
                                {
                                    rv[lr].emplace_back(tree.left, tree.right);
                                    lr = tree.right_sib[lr];
                                }
                        }
                }
            return rv;
        }
    } // namespace ts
} // namespace fwdpp

#endif
