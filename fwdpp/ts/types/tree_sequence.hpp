#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#if __cplusplus >= 201703L
#include <optional>
#endif
#include "table_collection.hpp"
#include "../marginal_tree.hpp"
#include "../exceptions.hpp"
#include "../node_flags.hpp"
#include "../tree_flags.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            template <typename SignedInteger> class tree_sequence;

            template <typename SignedInteger> class tree_iterator
            {
              private:
                std::reference_wrapper<const tree_sequence<SignedInteger>> treeseq_;
                marginal_tree<SignedInteger> tree_;
                std::size_t input_edge_index, output_edge_index;
                double position;

                bool
                still_advancing() const
                // assumes left-to-right iteration
                {
                    return (input_edge_index < tables().input_left.size()
                            || position < tables().genome_length());
                }

              public:
                tree_iterator(const tree_sequence<SignedInteger>& ts,
                              std::uint32_t flags)
                    : treeseq_{ts}, tree_{ts.num_nodes(), ts.samples(),
                                          static_cast<bool>(
                                              flags & tree_flags::TRACK_SAMPLES)},
                      input_edge_index{0}, output_edge_index{0}, position{0.}
                {
                }

                const marginal_tree<SignedInteger>*
                tree_ptr() const
                {
                    return still_advancing() ? &tree_ : nullptr;
                }

#if __cplusplus >= 201703L
                std::optional<std::reference_wrapper<const marginal_tree<SignedInteger>>>
                tree() const
                {
                    if (still_advancing())
                        {
                            return std::optional<std::reference_wrapper<
                                const marginal_tree<SignedInteger>>>{std::cref(tree_)};
                        }

                    return std::nullopt;
                }
#endif
                const types::table_collection<SignedInteger>&
                tables() const
                {
                    return treeseq_.get().tables();
                }
            };

            template <typename SignedInteger> class tree_sequence
            {
              private:
                std::shared_ptr<const types::table_collection<SignedInteger>> tables_;
                std::vector<SignedInteger> samples_;
                std::size_t num_trees_;

                static constexpr SignedInteger null
                    = types::table_collection<SignedInteger>::null;

                std::shared_ptr<const types::table_collection<SignedInteger>>
                init_tables(
                    std::shared_ptr<const types::table_collection<SignedInteger>> tables)
                {
                    if (tables == nullptr)
                        {
                            throw std::invalid_argument(
                                "input pointer to table_collection is nullptr");
                        }
                    if (!tables->indexed())
                        {
                            throw tables_error("tables are not indexed");
                        }
                    return tables;
                }

                std::vector<SignedInteger>
                init_samples(const types::table_collection<SignedInteger>& tables)
                {
                    std::vector<SignedInteger> rv;
                    for (std::size_t i = 0; i < tables.nodes.size(); ++i)
                        {
                            if (tables.nodes[i].flags & node_flags::IS_SAMPLE)
                                {
                                    rv.push_back(static_cast<SignedInteger>(i));
                                }
                        }
                    return rv;
                }

                std::vector<SignedInteger>
                init_samples(const types::table_collection<SignedInteger>& tables,
                             std::vector<SignedInteger> samples)
                {
                    std::vector<std::size_t> is_sample(tables.nodes.size(), 0);
                    for (auto s : samples)
                        {
                            if (s == null)
                                {
                                    throw std::invalid_argument(
                                        "input sample ID is null");
                                }
                            if (is_sample[static_cast<std::size_t>(s)] == 1)
                                {
                                    throw std::invalid_argument(
                                        "node ID marked as sample multiple times");
                                }
                            is_sample[static_cast<std::size_t>(s)] = 1;
                        }
                    return samples;
                }

                std::size_t
                init_num_trees()
                {
                    return 0;
                }

              public:
                tree_sequence(
                    std::shared_ptr<types::table_collection<SignedInteger>> tables)
                    : tables_{init_tables(std::move(tables))},
                      samples_{init_samples(*tables)}, num_trees_{init_num_trees()}
                {
                }

                tree_sequence(
                    std::shared_ptr<types::table_collection<SignedInteger>> tables,
                    std::vector<SignedInteger> samples)
                    : tables_{init_tables(std::move(tables))}, samples_{init_samples(
                                                                   *tables,
                                                                   std::move(samples))},
                      num_trees_{init_num_trees()}
                {
                }

                std::size_t
                num_nodes() const
                {
                    return tables_->nodes.size();
                }

                const std::vector<SignedInteger>&
                samples() const
                {
                    return samples_;
                }

                tree_iterator<SignedInteger>
                trees(std::uint32_t flags)
                {
                    return tree_iterator<SignedInteger>{*this, flags};
                }

                const types::table_collection<SignedInteger>&
                tables() const
                {
                    return *tables_;
                }
            };

#if __cplusplus < 201703L
            template <typename SignedInteger>
            constexpr SignedInteger tree_sequence<SignedInteger>::null;
#endif
        }
    }
}
