#ifndef FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP
#define FWDPP_TS_MARK_MULTIPLE_ROOTS_HPP

#include <map>
#include <utility>
#include <vector>
#include "definitions.hpp"
#include "table_collection.hpp"
#include "marginal_tree_iterator.hpp"

namespace fwdpp
{
    namespace ts
    {
        // TODO: consider flattening the return value to a vector
        std::map<TS_NODE_INT, std::vector<std::pair<double, double>>>
        mark_multiple_roots(const table_collection &tables,
                            const std::vector<TS_NODE_INT> &samples)
        {
            std::map<TS_NODE_INT, std::vector<std::pair<double, double>>> rv;
            marginal_tree_iterator mti(tables, samples);
            while (mti(std::true_type(), std::false_type()))
                {
                    bool single_root = false;
					auto & tree = mti.tree();
                    for (auto &s : samples)
                        {
                            auto p = s;
                            auto lp = p;
                            while (p != -1)
                                {
                                    lp = p;
                                    p = tree.parents[p];
                                }
                            if (tree.leaf_counts[lp] == samples.size())
                                {
                                    single_root = true;
                                }
                            else
                                {
                                    auto itr = rv.find(lp);
                                    auto w = std::make_pair(
                                        tree.left, tree.right);
                                    if (itr == rv.end())
                                        {
                                            rv[lp].emplace_back(std::move(w));
                                        }
                                    else
                                        {
                                            assert(!itr->second.empty());
                                            if (std::find(itr->second.begin(),
                                                          itr->second.end(), w)
                                                == itr->second.end())
                                                {
                                                    itr->second.emplace_back(
                                                        std::move(w));
                                                }
                                        }
                                }
                            if (single_root)
                                {
                                    break;
                                }
                        }
                }
            return rv;
        }
    } // namespace ts
} // namespace fwdpp

#endif
