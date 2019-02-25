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
        /// \version 0.7.4 Update to include data structures for root tracking
        {
          private:
            void
            init_samples(const std::vector<TS_NODE_INT>& samples)
            {
                for (std::size_t i = 0; i < samples.size(); ++i)
                    {
                        auto s = samples[i];
                        // See GitHub issue #158 for background
                        if (sample_index_map[s] != TS_NULL_NODE)
                            {
                                throw std::invalid_argument(
                                    "invalid sample list");
                            }
                        sample_index_map[s] = i;
                        ++sample_size;
                        left_sample[s] = right_sample[s] = sample_index_map[s];
                        above_sample[s] = 1;
                        // Initialize roots
                        if (i < samples.size() - 1)
                            {
                                right_sib[s] = samples[i + 1];
                            }
                        if (i > 0)
                            {
                                left_sib[s] = samples[i - 1];
                            }
                    }
            }

          public:
            std::vector<TS_NODE_INT> parents, leaf_counts,
                preserved_leaf_counts, left_sib, right_sib, left_child,
                right_child, left_sample, right_sample, next_sample,
                sample_index_map;
            std::vector<std::int8_t> above_sample;
            double left, right;
            TS_NODE_INT sample_size, left_root;
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
                  above_sample(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() },
                  right{ std::numeric_limits<double>::quiet_NaN() },
                  sample_size(0), left_root(TS_NULL_NODE)
            /// Constructor
            /// \todo Document
            {
                init_samples(samples);
                for (auto s : samples)
                    {
                        leaf_counts[s] = 1;
                    }
                left_root = samples[0];
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
                  above_sample(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() },
                  right{ std::numeric_limits<double>::quiet_NaN() },
                  sample_size(0), left_root(TS_NULL_NODE)
            /// Constructor
            /// \todo Document
            {
                auto merged_samples(samples);
                if (!preserved_nodes.empty())
                    {
                        merged_samples.insert(merged_samples.end(),
                                              begin(preserved_nodes),
                                              end(preserved_nodes));
                    }
                init_samples(merged_samples);
                left_root = samples[0];
                for (auto s : samples)
                    {
                        leaf_counts[s] = 1;
                    }
                for (auto s : preserved_nodes)
                    {
                        preserved_leaf_counts[s] = 1;
                    }
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
                  above_sample(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() },
                  right{ std::numeric_limits<double>::quiet_NaN() },
                  sample_size(0), left_root(TS_NULL_NODE)
            /// Constructor
            /// \todo Document
            {
            }

            int
            num_roots() const
            /// Return number of roots
            {
                if (left_root == TS_NULL_NODE)
                    {
                        throw std::runtime_error("left_root is NULL");
                    }
                int nroots = 0;
                auto lr = left_root;
                while (lr != TS_NULL_NODE)
                    {
                        ++nroots;
                        lr = right_sib[lr];
                    }
                return nroots;
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
