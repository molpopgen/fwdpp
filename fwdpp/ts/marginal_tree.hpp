#ifndef FWDPP_TS_MARGINAL_TREE_HPP
#define FWDPP_TS_MARGINAL_TREE_HPP

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <limits>
#include <cstdint>
#include "definitions.hpp"
#include "exceptions.hpp"

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
        /// \version 0.8.0 Now holds a list of samples.
        {
          private:
            std::vector<TS_NODE_INT> samples_list;

            std::vector<TS_NODE_INT>
            forward_input_samples(std::vector<TS_NODE_INT>&& a,
                                  std::vector<TS_NODE_INT>&& b)
            {
                std::vector<TS_NODE_INT> rv(std::move(a));
                std::move(begin(b), end(b), std::back_inserter(a));
                return rv;
            }

            std::vector<TS_NODE_INT>
            forward_input_samples(const std::vector<TS_NODE_INT>& a,
                                  const std::vector<TS_NODE_INT>& b)
            {
                std::vector<TS_NODE_INT> rv(a);
                rv.insert(rv.end(), begin(b), end(b));
                return rv;
            }

            void
            init_samples()
            {
                for (std::size_t i = 0; i < samples_list.size(); ++i)
                    {
                        auto s = samples_list[i];
                        // See GitHub issue #158 for background
                        if (sample_index_map[s] != TS_NULL_NODE)
                            {
                                throw std::invalid_argument(
                                    "invalid sample list");
                            }
                        sample_index_map[s] = i;
                        left_sample[s] = right_sample[s] = sample_index_map[s];
                        above_sample[s] = 1;
                        // Initialize roots
                        if (i < samples_list.size() - 1)
                            {
                                right_sib[s] = samples_list[i + 1];
                            }
                        if (i > 0)
                            {
                                left_sib[s] = samples_list[i - 1];
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
            TS_NODE_INT left_root;

            template <typename SAMPLES>
            marginal_tree(TS_NODE_INT nnodes, SAMPLES&& samples)
                : samples_list(std::forward<SAMPLES>(samples)),
                  parents(nnodes, TS_NULL_NODE), leaf_counts(nnodes, 0),
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
                  left_root(TS_NULL_NODE)
            /// Constructor
            /// \todo Document
            {
                if (samples_list.empty())
                    {
                        throw empty_samples(
                            "marginal_tree: empty sample list");
                    }
                init_samples();
                for (auto s : samples_list)
                    {
                        leaf_counts[s] = 1;
                    }
                left_root = samples_list[0];
            }

            template <typename SAMPLES>
            marginal_tree(TS_NODE_INT nnodes, SAMPLES&& samples,
                          SAMPLES&& preserved_nodes)
                : samples_list(forward_input_samples(
                    std::forward<SAMPLES>(samples),
                    std::forward<SAMPLES>(preserved_nodes))),
                  parents(nnodes, TS_NULL_NODE), leaf_counts(nnodes, 0),
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
                  left_root(TS_NULL_NODE)
            /// Constructor
            /// \todo Document
            {
                if (samples_list.empty())
                    {
                        throw empty_samples(
                            "marginal_tree: empty sample list");
                    }
                init_samples();
                left_root = samples_list[0];
                for (auto s : samples_list)
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
                  left_root(TS_NULL_NODE)
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

            inline std::size_t
            sample_size() const
            {
                return samples_list.size();
            }

            inline std::vector<TS_NODE_INT>::const_iterator
            samples_list_begin() const
            {
                return begin(samples_list);
            }

            inline std::vector<TS_NODE_INT>::const_iterator
            samples_list_end() const
            {
                return end(samples_list);
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
