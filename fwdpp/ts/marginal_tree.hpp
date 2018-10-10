#ifndef FWDPP_TS_MARGINAL_TREE_HPP
#define FWDPP_TS_MARGINAL_TREE_HPP

#include <stdexcept>
#include <vector>
#include <limits>
#include <cstdint>
#include "definitions.hpp"

namespace fwdpp
{
    namespace ts
    {
        struct marginal_tree
        {
            std::vector<std::int32_t> parents, leaf_counts,
                preserved_leaf_counts, left_sib, right_sib, left_child,
                right_child, left_sample, right_sample, next_sample, is_sample;
            double left, right;
            marginal_tree(std::int32_t nnodes,
                          const std::vector<std::int32_t>& samples)
                : parents(nnodes, TS_NULL_NODE), leaf_counts(nnodes, 0),
                  preserved_leaf_counts(nnodes, 0), left_sib(nnodes, TS_NULL_NODE),
                  right_sib(nnodes, TS_NULL_NODE), left_child(nnodes, TS_NULL_NODE),
                  right_child(nnodes, TS_NULL_NODE), left_sample(nnodes, TS_NULL_NODE),
                  right_sample(nnodes, TS_NULL_NODE), next_sample(nnodes, TS_NULL_NODE),
                  is_sample(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  }
            {
                for (auto s : samples)
                    {
                        if (static_cast<std::size_t>(s) >= leaf_counts.size())
                            {
                                throw std::invalid_argument(
                                    "sample index out of range");
                            }
                        leaf_counts[s] = 1;
                        is_sample[s] = 1;
                        left_sample[s] = right_sample[s] = s;
                    }
            }
            marginal_tree(std::int32_t nnodes,
                          const std::vector<std::int32_t>& samples,
                          const std::vector<std::int32_t>& preserved_nodes)
                : parents(nnodes, TS_NULL_NODE), leaf_counts(nnodes, 0),
                  preserved_leaf_counts(nnodes, 0), left_sib(nnodes, TS_NULL_NODE),
                  right_sib(nnodes, TS_NULL_NODE), left_child(nnodes, TS_NULL_NODE),
                  right_child(nnodes, TS_NULL_NODE), left_sample(nnodes, TS_NULL_NODE),
                  right_sample(nnodes, TS_NULL_NODE), next_sample(nnodes, TS_NULL_NODE),
                  is_sample(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  }
            {
                for (auto s : samples)
                    {
                        if (static_cast<std::size_t>(s) >= leaf_counts.size())
                            {
                                throw std::invalid_argument(
                                    "sample index out of range");
                            }
                        leaf_counts[s] = 1;
                        is_sample[s] = 1;
                        left_sample[s] = right_sample[s] = s;
                    }
                for (auto s : preserved_nodes)
                    {
                        if (static_cast<std::size_t>(s) >= leaf_counts.size())
                            {
                                throw std::invalid_argument(
                                    "sample index out of range");
                            }
                        preserved_leaf_counts[s] = 1;
                        is_sample[s] = 1;
                        left_sample[s] = right_sample[s] = s;
                    }
            }
            marginal_tree(std::int32_t nnodes)
                : parents(nnodes, TS_NULL_NODE), leaf_counts{}, preserved_leaf_counts{},
                  left_sib(nnodes, TS_NULL_NODE), right_sib(nnodes, TS_NULL_NODE),
                  left_child(nnodes, TS_NULL_NODE), right_child(nnodes, TS_NULL_NODE),
                  left_sample(nnodes, TS_NULL_NODE), right_sample(nnodes, TS_NULL_NODE),
                  next_sample(nnodes, TS_NULL_NODE), is_sample(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  }
            {
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
