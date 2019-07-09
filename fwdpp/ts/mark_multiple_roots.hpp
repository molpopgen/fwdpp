#ifndef FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP
#define FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP

#include <map>
#include <utility>
#include <vector>
#include "definitions.hpp"
#include "table_collection.hpp"
#include "tree_visitor.hpp"
#include "marginal_tree_functions/roots.hpp"

namespace fwdpp
{
    namespace ts
    {
        // TODO: consider flattening the return value to a vector
        template <typename SAMPLES>
        inline std::map<TS_NODE_INT, std::vector<std::pair<double, double>>>
        mark_multiple_roots(const table_collection &tables, SAMPLES &&samples)
        /// \brief Identify root nodes in "marginal forests".
        ///
        /// \version 0.7.0 Added to library
        /// \version 0.7.4 Refactored to use root tracking method
        /// \version 0.8.0 Refactored to use ts::num_roots and ts::process_roots
        ///
        /// See fwdpp::ts::mutate_tables for discussion.
        {
            std::map<TS_NODE_INT, std::vector<std::pair<double, double>>> rv;
            tree_visitor mti(tables, std::forward<SAMPLES>(samples),
                             update_samples_list(false));
            while (mti())
                {
                    auto &tree = mti.tree();
                    if (num_roots(tree) > 1)
                        {
                            const auto func = [&rv, &tree](TS_NODE_INT root) {
                                rv[root].emplace_back(tree.left, tree.right);
                            };
                            process_roots(tree, func);
                        }
                }
            return rv;
        }
    } // namespace ts
} // namespace fwdpp

#endif
