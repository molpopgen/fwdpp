#ifndef FWDPP_TS_MARGINAL_TREE_HPP
#define FWDPP_TS_MARGINAL_TREE_HPP

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <limits>
#include <cstdint>
#include <type_traits>
#include "exceptions.hpp"
#include "types/generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger = std::int32_t> struct sample_group_map
        /// \brief Maps a node id to a sample group
        ///
        /// When constructing a fwdpp::ts::tree_visitor,
        /// vectors of this type may be used to mark
        /// sample nodes as beloning to different "groups".
        {
            static_assert(std::is_integral<SignedInteger>::value,
                          "SignedInteger must be an integral type");
            static_assert(std::is_signed<SignedInteger>::value,
                          "SignedInteger must be a signed type");
            /// node id
            SignedInteger node_id;
            /// group label
            std::int32_t group;
            sample_group_map(SignedInteger n, std::int32_t g) : node_id(n), group(g)
            /// \param n A node id
            /// \param g A gropup label
            {
                if (n == types::generate_null_id<SignedInteger>())
                    {
                        throw samples_error("null ID passed to sample_index_map");
                    }
            }
        };

        template <typename SignedInteger = std::int32_t> class marginal_tree
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
        /// \version 0.8.0 Now holds a list of samples. Samples may be assigned to groups.
        {
          private:
            static_assert(std::is_integral<SignedInteger>::value,
                          "SignedInteger must be an integral type");
            static_assert(std::is_signed<SignedInteger>::value,
                          "SignedInteger must be a signed type");
            std::size_t num_nodes;
            std::vector<std::int32_t> sample_groups;
            std::vector<SignedInteger> samples_list;
            bool advancing_sample_list_;

            static constexpr SignedInteger null
                = types::generate_null_id<SignedInteger>();

            std::vector<std::int32_t>
            fill_sample_groups(const std::vector<SignedInteger>& samples)
            {
                std::vector<std::int32_t> rv(num_nodes,
                                             std::numeric_limits<std::int32_t>::min());
                for (auto i : samples)
                    {
                        rv[i] = 0;
                    }
                return rv;
            }

            std::vector<std::int32_t>
            fill_sample_groups(
                const std::vector<sample_group_map<SignedInteger>>& samples)
            {
                std::vector<std::int32_t> rv(num_nodes,
                                             std::numeric_limits<std::int32_t>::min());
                for (auto i : samples)
                    {
                        rv[i.node_id] = i.group;
                    }
                return rv;
            }

            std::vector<std::int32_t>
            fill_sample_groups(const std::vector<SignedInteger>& samples_a,
                               const std::vector<SignedInteger>& samples_b)
            {
                std::vector<std::int32_t> rv(num_nodes,
                                             std::numeric_limits<std::int32_t>::min());
                for (auto i : samples_a)
                    {
                        rv[i] = 0;
                    }
                for (auto i : samples_b)
                    {
                        rv[i] = 1;
                    }
                return rv;
            }

            std::vector<std::int32_t>
            fill_sample_groups(
                const std::vector<sample_group_map<SignedInteger>>& samples_a,
                const std::vector<SignedInteger>& samples_b)
            {
                std::vector<std::int32_t> rv(num_nodes,
                                             std::numeric_limits<std::int32_t>::min());
                for (auto i : samples_a)
                    {
                        rv[i.node_id] = 0;
                    }
                for (auto i : samples_b)
                    {
                        rv[i] = 1;
                    }
                return rv;
            }

            std::vector<SignedInteger>
            init_samples_list(const std::vector<SignedInteger>& s)
            {
                return s;
            }

            std::vector<SignedInteger>
            init_samples_list(const std::vector<sample_group_map<SignedInteger>>& s)
            {
                std::vector<SignedInteger> rv;
                for (auto& i : s)
                    {
                        rv.push_back(i.node_id);
                    }
                return rv;
            }

            std::vector<SignedInteger>
            init_samples_list(const std::vector<SignedInteger>& a,
                              const std::vector<SignedInteger>& b)
            {
                auto rv = a;
                rv.insert(end(rv), begin(b), end(b));
                return rv;
            }

            std::vector<SignedInteger>
            init_samples_list(const std::vector<sample_group_map<SignedInteger>>& a,
                              const std::vector<SignedInteger>& b)
            {
                auto rv = init_samples_list(a);
                rv.insert(end(rv), begin(b), end(b));
                return rv;
            }

            void
            init_samples()
            {
                for (std::size_t i = 0; i < samples_list.size(); ++i)
                    {
                        auto s = samples_list[i];
                        // See GitHub issue #158 for background
                        if (sample_index_map[s] != null)
                            {
                                throw samples_error("invalid sample list");
                            }
                        sample_index_map[s] = static_cast<SignedInteger>(i);
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
            using id_type = SignedInteger;
            std::vector<id_type> parents, leaf_counts, preserved_leaf_counts, left_sib,
                right_sib, left_child, right_child, left_sample, right_sample,
                next_sample, sample_index_map;
            std::vector<std::int8_t> above_sample;
            double left, right;
            id_type left_root;

            template <typename SAMPLES>
            marginal_tree(std::size_t nnodes, const SAMPLES& samples,
                          bool advancing_sample_list)
                : num_nodes(nnodes), sample_groups(fill_sample_groups(samples)),
                  samples_list(init_samples_list(samples)),
                  advancing_sample_list_(advancing_sample_list), parents(nnodes, null),
                  leaf_counts(nnodes, 0), preserved_leaf_counts(nnodes, 0),
                  left_sib(nnodes, null), right_sib(nnodes, null),
                  left_child(nnodes, null), right_child(nnodes, null),
                  left_sample(nnodes, null), right_sample(nnodes, null),
                  next_sample(nnodes, null), sample_index_map(nnodes, null),
                  above_sample(nnodes, 0),
                  left{std::numeric_limits<double>::quiet_NaN()},
                  right{std::numeric_limits<double>::quiet_NaN()}, left_root(null)
            /// \param nnodes Number of nodes in table_collection
            /// \param samples The sample list
            ///
            /// \a samples may be either std::vector<id_type> or
            /// std::vector<fwdpp::ts::sample_group_map>.  For the former, all
            /// sample nodes will be assigned group 0.
            {
                if (samples_list.empty())
                    {
                        throw samples_error("marginal_tree: empty sample list");
                    }
                init_samples();
                for (auto s : samples_list)
                    {
                        leaf_counts[s] = 1;
                    }
                left_root = samples_list[0];
            }

            template <typename SAMPLES>
            marginal_tree(std::size_t nnodes, const SAMPLES& samples,
                          const std::vector<id_type> preserved_nodes,
                          bool advancing_sample_list)
                : num_nodes(nnodes),
                  sample_groups(fill_sample_groups(samples, preserved_nodes)),
                  samples_list(init_samples_list(samples, preserved_nodes)),
                  advancing_sample_list_(advancing_sample_list), parents(nnodes, null),
                  leaf_counts(nnodes, 0), preserved_leaf_counts(nnodes, 0),
                  left_sib(nnodes, null), right_sib(nnodes, null),
                  left_child(nnodes, null), right_child(nnodes, null),
                  left_sample(nnodes, null), right_sample(nnodes, null),
                  next_sample(nnodes, null), sample_index_map(nnodes, null),
                  above_sample(nnodes, 0),
                  left{std::numeric_limits<double>::quiet_NaN()},
                  right{std::numeric_limits<double>::quiet_NaN()}, left_root(null)
            /// Constructor
            {
                if (samples_list.empty())
                    {
                        throw samples_error("marginal_tree: empty sample list");
                    }
                init_samples();
                left_root = samples_list[0];
                for (std::size_t i = 0; i < samples.size(); ++i)
                    {
                        leaf_counts[samples_list[i]] = 1;
                    }
                for (auto s : preserved_nodes)
                    {
                        preserved_leaf_counts[s] = 1;
                    }
            }

            marginal_tree(id_type nnodes)
                : num_nodes(nnodes), sample_groups{}, samples_list{},
                  advancing_sample_list_(false),
                  parents(nnodes, null), leaf_counts{}, preserved_leaf_counts{},
                  left_sib(nnodes, null), right_sib(nnodes, null),
                  left_child(nnodes, null), right_child(nnodes, null),
                  left_sample(nnodes, null), right_sample(nnodes, null),
                  next_sample(nnodes, null), sample_index_map(nnodes, null),
                  above_sample(nnodes, 0),
                  left{std::numeric_limits<double>::quiet_NaN()},
                  right{std::numeric_limits<double>::quiet_NaN()}, left_root(null)
            /// Constructor
            /// \todo Document
            {
            }

            int
            num_roots() const
            /// Return number of roots
            {
                if (left_root == null)
                    {
                        throw std::runtime_error("left_root is NULL");
                    }
                int nroots = 0;
                auto lr = left_root;
                while (lr != null)
                    {
                        ++nroots;
                        lr = right_sib[lr];
                    }
                return nroots;
            }

            inline std::size_t
            sample_size() const
            /// Number of samples
            {
                return samples_list.size();
            }

            inline typename std::vector<id_type>::const_iterator
            samples_list_begin() const
            /// Beginning of samples list
            {
                return begin(samples_list);
            }

            inline typename std::vector<id_type>::const_iterator
            samples_list_end() const
            /// End itertor for samples list
            {
                return end(samples_list);
            }

            inline id_type
            sample_group(id_type u) const
            /// Return the sample group for node \u
            {
                if (static_cast<std::size_t>(u) >= num_nodes)
                    {
                        throw std::invalid_argument("invalid node id");
                    }
                return sample_groups[u];
            }

            inline bool
            advancing_sample_list() const
            {
                return advancing_sample_list_;
            }

            inline std::size_t
            size() const
            /// Return the length of the internal vectors.
            {
                return num_nodes;
            }

            inline id_type
            sample_table_index_to_node(id_type u) const
            /// If u is a sample index, return the associated
            /// node id.
            {
                return samples_list[u];
            }
        };

#if __cplusplus < 201703L
            template <typename SignedInteger>
            constexpr SignedInteger marginal_tree<SignedInteger>::null;
#endif

    } // namespace ts
} // namespace fwdpp

#endif
