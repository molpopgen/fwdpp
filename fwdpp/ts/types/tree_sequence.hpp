#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include "table_collection.hpp"
#include "../marginal_tree.hpp"
#include "../exceptions.hpp"
#include "../node_flags.hpp"
#include "../tree_sequence_flags.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            template <typename SignedInteger> class tree_sequence
            {
              private:
                std::shared_ptr<const types::table_collection<SignedInteger>> tables_;
                std::vector<SignedInteger> samples_;
                std::size_t num_trees_;
                std::uint32_t flags_;
                marginal_tree<SignedInteger> tree_;

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
                    std::shared_ptr<types::table_collection<SignedInteger>> tables,
                    std::uint32_t flags)
                    : tables_{init_tables(std::move(tables))}, samples_{init_samples(
                                                                   *tables)},
                      num_trees_{init_num_trees()}, flags_{flags},
                      tree_{
                          tables->nodes.size(), samples_,
                          static_cast<bool>(flags_ & tree_sequence_flags::TRACK_SAMPLES)}
                {
                }

                tree_sequence(
                    std::shared_ptr<types::table_collection<SignedInteger>> tables,
                    std::vector<SignedInteger> samples, std::uint32_t flags)
                    : tables_{init_tables(std::move(tables))}, samples_{init_samples(
                                                                   *tables,
                                                                   std::move(samples))},
                      num_trees_{init_num_trees()}, flags_{flags},
                      tree_{
                          tables->nodes.size(), samples_,
                          static_cast<bool>(flags_ & tree_sequence_flags::TRACK_SAMPLES)}
                {
                }
            };

#if __cplusplus < 201703L
            template <typename SignedInteger>
            constexpr SignedInteger tree_sequence<SignedInteger>::null;
#endif
        }
    }
}
