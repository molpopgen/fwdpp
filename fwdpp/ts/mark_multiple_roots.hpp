#ifndef FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP
#define FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP

#include <map>
#include <utility>
#include <vector>
#include "definitions.hpp"
#include "tree_visitor.hpp"
#include "marginal_tree_functions/roots.hpp"

namespace fwdpp
{
    namespace ts
    {
        // TODO: consider flattening the return value to a vector
        template <typename TableCollectionType, typename SAMPLES>
        inline std::map<table_index_t, std::vector<std::pair<double, double>>>
        mark_multiple_roots(const TableCollectionType &tables, SAMPLES &&samples)
        /// \brief Identify root nodes in "marginal forests".
        ///
        /// \version 0.7.0 Added to library
        /// \version 0.7.4 Refactored to use root tracking method
        /// \version 0.8.0 Refactored to use ts::num_roots and ts::process_roots
        /// \version 0.9.0 Added typename TableCollectionType
        ///
        /// See fwdpp::ts::mutate_tables for discussion.
        {
            std::map<table_index_t, std::vector<std::pair<double, double>>> rv;
            tree_visitor<TableCollectionType> mti(tables, std::forward<SAMPLES>(samples),
                                                  update_samples_list(false));
            while (mti())
                {
                    auto &tree = mti.tree();
                    if (num_roots(tree) > 1)
                        {
                            const auto func = [&rv, &tree](table_index_t root) {
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
