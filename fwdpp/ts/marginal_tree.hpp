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
        class marginal_tree
        /// \brief A non-recombining tree
        ///
        /// The tree is represented as a sparse tree
        /// data structure where several linked lists
        /// allow efficient retrieval of parent/child/siblings
        /// of nodes.
        ///
        /// \version 0.7.0 Added to fwdpp
        /// \version 0.7.1 Constructors throw exceptions when sample lists contain the
        /// same node ID more than once. Changed from struct to class in order
        /// to reuse some code in a private function. Initialization also 
        /// tracks the total sample size, which is number of nonzero elements
        /// in sample_index_map.
        {
          private:
            void
            init_samples(const std::vector<TS_NODE_INT>& samples,
                         std::vector<TS_NODE_INT>& lc)
            {
                int i = 0;
                for (auto s : samples)
                    {
                        if (static_cast<std::size_t>(s) >= leaf_counts.size())
                            {
                                throw std::invalid_argument(
                                    "sample index out of range");
                            }
                        // See GitHub issue #158 for background
                        if (sample_index_map[s] != TS_NULL_NODE)
                            {
                                throw std::invalid_argument(
                                    "invalid sample list");
                            }
                        lc[s] = 1;
                        sample_index_map[s] = i;
                        ++sample_size;
                        left_sample[s] = right_sample[s] = sample_index_map[s];
                        i++;
                    }
            }

          public:
            std::vector<TS_NODE_INT> parents, leaf_counts,
                preserved_leaf_counts, left_sib, right_sib, left_child,
                right_child, left_sample, right_sample, next_sample,
                sample_index_map;
            double left, right;
            TS_NODE_INT sample_size;
            marginal_tree(TS_NODE_INT nnodes,
                          const std::vector<TS_NODE_INT>& samples)
                : parents(nnodes, TS_NULL_NODE), leaf_counts(nnodes, 0),
                  preserved_leaf_counts(nnodes, 0),
                  left_sib(nnodes, TS_NULL_NODE),
                  right_sib(nnodes, TS_NULL_NODE),
                  left_child(nnodes, TS_NULL_NODE),
                  right_child(nnodes, TS_NULL_NODE),
                  left_sample(nnodes, TS_NULL_NODE),
                  right_sample(nnodes, TS_NULL_NODE),
                  next_sample(nnodes, TS_NULL_NODE),
                  sample_index_map(nnodes, TS_NULL_NODE),
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  },sample_size(0)
            /// Constructor
            /// \todo Document
            {
                init_samples(samples, leaf_counts);
            }

            marginal_tree(TS_NODE_INT nnodes,
                          const std::vector<TS_NODE_INT>& samples,
                          const std::vector<TS_NODE_INT>& preserved_nodes)
                : parents(nnodes, TS_NULL_NODE), leaf_counts(nnodes, 0),
                  preserved_leaf_counts(nnodes, 0),
                  left_sib(nnodes, TS_NULL_NODE),
                  right_sib(nnodes, TS_NULL_NODE),
                  left_child(nnodes, TS_NULL_NODE),
                  right_child(nnodes, TS_NULL_NODE),
                  left_sample(nnodes, TS_NULL_NODE),
                  right_sample(nnodes, TS_NULL_NODE),
                  next_sample(nnodes, TS_NULL_NODE),
                  sample_index_map(nnodes, TS_NULL_NODE),
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  },sample_size(0)
            /// Constructor
            /// \todo Document
            {
                init_samples(samples, leaf_counts);
                init_samples(preserved_nodes, preserved_leaf_counts);
            }
            marginal_tree(TS_NODE_INT nnodes)
                : parents(nnodes, TS_NULL_NODE), leaf_counts{},
                  preserved_leaf_counts{}, left_sib(nnodes, TS_NULL_NODE),
                  right_sib(nnodes, TS_NULL_NODE),
                  left_child(nnodes, TS_NULL_NODE),
                  right_child(nnodes, TS_NULL_NODE),
                  left_sample(nnodes, TS_NULL_NODE),
                  right_sample(nnodes, TS_NULL_NODE),
                  next_sample(nnodes, TS_NULL_NODE),
                  sample_index_map(nnodes, TS_NULL_NODE),
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  },sample_size(0)
            /// Constructor
            /// \todo Document
            {
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
