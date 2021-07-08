#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <iterator>
#if __cplusplus >= 201703L
#include <optional>
#endif
#include "table_collection.hpp"
#include "../marginal_tree.hpp"
#include "../exceptions.hpp"
#include "../node_flags.hpp"
#include "../tree_flags.hpp"
#include "../detail/advance_marginal_tree_policies.hpp"

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
                // NOTE: since tree stores the samples,
                // we can probably just keep a shared_ptr to the tables here,
                // adding a level of memory safety.
                std::reference_wrapper<const tree_sequence<SignedInteger>> treeseq_;
                marginal_tree<SignedInteger> tree_;
                double position, maxpos;
                bool advanced_;
                typename std::vector<SignedInteger>::const_iterator current_input_edge,
                    input_edge_end, current_output_edge, output_edge_end;
                typename types::table_collection<
                    SignedInteger>::edge_table::const_iterator beg_edges,
                    end_edges;

              public:
                tree_iterator(const tree_sequence<SignedInteger>& ts,
                              std::uint32_t flags)
                    : treeseq_{ts}, tree_{ts.num_nodes(), ts.samples(),
                                          static_cast<bool>(
                                              flags & tree_flags::TRACK_SAMPLES)},
                      position{0.}, maxpos{tables().genome_length()}, advanced_{false},
                      current_input_edge{end(tables().input_left)},
                      input_edge_end{end(tables().input_left)},
                      current_output_edge{end(tables().output_right)},
                      output_edge_end{end(tables().output_right)},
                      beg_edges{begin(tables().edges)}, end_edges{end(tables().edges)}
                {
                }

                const marginal_tree<SignedInteger>*
                tree_ptr() const
                // tree() should be preferred for c++ >= 17
                {
                    return advanced_ ? &tree_ : nullptr;
                }

#if __cplusplus >= 201703L
                std::optional<std::reference_wrapper<const marginal_tree<SignedInteger>>>
                tree() const
                {
                    if (advanced_)
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

                void
                advance_right()
                // Move 1 tree left-to-right
                {
                    if (current_input_edge < input_edge_end || position < maxpos)
                        {
                            while (current_output_edge < output_edge_end

                                   && (beg_edges + *current_output_edge)->right
                                          == position) // T4
                                {
                                    const auto p
                                        = (beg_edges + *current_output_edge)->parent;
                                    const auto c
                                        = (beg_edges + *current_output_edge)->child;
                                    const auto lsib = tree_.left_sib[c];
                                    const auto rsib = tree_.right_sib[c];
                                    if (lsib
                                        == types::table_collection<SignedInteger>::null)
                                        {
                                            tree_.left_child[p] = rsib;
                                        }
                                    else
                                        {
                                            tree_.right_sib[lsib] = rsib;
                                        }
                                    if (rsib
                                        == types::table_collection<SignedInteger>::null)
                                        {
                                            tree_.right_child[p] = lsib;
                                        }
                                    else
                                        {
                                            tree_.left_sib[rsib] = lsib;
                                        }
                                    tree_.parents[c]
                                        = types::table_collection<SignedInteger>::null;
                                    tree_.left_sib[c]
                                        = types::table_collection<SignedInteger>::null;
                                    tree_.right_sib[c]
                                        = types::table_collection<SignedInteger>::null;
                                    detail::outgoing_leaf_counts(
                                        tree_,
                                        (beg_edges + *current_output_edge)->parent,
                                        (beg_edges + *current_output_edge)->child);
                                    if (tree_.advancing_sample_list())
                                        {
                                            detail::update_samples_list(
                                                tree_, (beg_edges + *current_output_edge)
                                                           ->parent);
                                        }
                                    detail::update_roots_outgoing(p, c, tree_);
                                    ++current_input_edge;
                                }
                            while (current_input_edge < input_edge_end
                                   && (beg_edges + *current_input_edge)->left
                                          == position) // Step T2
                                {
                                    const auto p
                                        = (beg_edges + *current_input_edge)->parent;
                                    const auto c
                                        = (beg_edges + *current_input_edge)->child;
                                    const auto rchild = tree_.right_child[p];
                                    const auto lsib = tree_.left_sib[c];
                                    const auto rsib = tree_.right_sib[c];
                                    if (rchild
                                        == types::table_collection<SignedInteger>::null)
                                        {
                                            tree_.left_child[p] = c;
                                            tree_.left_sib[c] = types::table_collection<
                                                SignedInteger>::null;
                                            tree_.right_sib[c] = types::table_collection<
                                                SignedInteger>::null;
                                        }
                                    else
                                        {
                                            tree_.right_sib[rchild] = c;
                                            tree_.left_sib[c] = rchild;
                                            tree_.right_sib[c] = types::table_collection<
                                                SignedInteger>::null;
                                        }
                                    // The entry for the child refers to
                                    // the parent's location in the node table.
                                    tree_.parents[c]
                                        = (beg_edges + *current_input_edge)->parent;
                                    tree_.right_child[p] = c;
                                    detail::incoming_leaf_counts(
                                        tree_, (beg_edges + *current_input_edge)->parent,
                                        (beg_edges + *current_input_edge)->child);
                                    if (tree_.advancing_sample_list())
                                        {
                                            detail::update_samples_list(
                                                tree_, (beg_edges + *current_input_edge)
                                                           ->parent);
                                        }
                                    detail::update_roots_incoming(p, c, lsib, rsib,
                                                                  tree_);

                                    ++current_input_edge;
                                }

                            // This is a big "gotcha".
                            // The root tracking functions will sometimes
                            // result in left_root not actually being the left_root.
                            // We loop through the left_sibs to fix that.
                            if (tree_.left_root
                                != types::table_collection<SignedInteger>::null)
                                {
                                    while (
                                        tree_.left_sib[tree_.left_root]
                                        != types::table_collection<SignedInteger>::null)
                                        {
                                            tree_.left_root
                                                = tree_.left_sib[tree_.left_root];
                                        }
                                }
#ifndef NDEBUG
                            // Validate the roots via brute-force.
                            auto lr = tree_.left_root;
                            if (lr == types::table_collection<SignedInteger>::null)
                                {
                                    throw std::runtime_error(
                                        "FWDPP DEBUG: left_root is null");
                                }
                            std::vector<int> is_root(tree_.sample_index_map.size(), 0);
                            std::vector<int> processed(is_root.size(), 0);
                            for (std::size_t s = 0; s < tree_.sample_index_map.size();
                                 ++s)
                                {
                                    if (tree_.sample_index_map[s]
                                        != types::table_collection<SignedInteger>::null)
                                        {
                                            auto u = static_cast<SignedInteger>(s);
                                            auto root = u;
                                            bool early_exit = false;
                                            while (u
                                                   != types::table_collection<
                                                       SignedInteger>::null)
                                                {
                                                    if (processed[u])
                                                        {
                                                            early_exit = true;
                                                            break;
                                                        }
                                                    processed[u] = 1;
                                                    root = u;
                                                    u = tree_.parents[u];
                                                }
                                            if (early_exit == false)
                                                {
                                                    is_root[root] = 1;
                                                }
                                        }
                                }
                            int nroots_brute = 0;
                            for (auto r : is_root)
                                {
                                    nroots_brute += r;
                                }
                            if (nroots_brute != tree_.num_roots())
                                {
                                    throw std::runtime_error("FWDPP DEBUG: num_roots "
                                                             "disagreement");
                                }
                            while (lr != types::table_collection<SignedInteger>::null)
                                {
                                    if (is_root[lr] != 1)
                                        {
                                            throw std::runtime_error("FWDPP DEBUG: root "
                                                                     "contents "
                                                                     "disagreement");
                                        }
                                    lr = tree_.right_sib[lr];
                                }
#endif
                            double right = maxpos;
                            if (current_input_edge < input_edge_end)
                                {
                                    right = std::min(
                                        right, (beg_edges + *current_input_edge)->left);
                                }
                            if (current_output_edge < output_edge_end)
                                {
                                    right = std::min(
                                        right,
                                        (beg_edges + *current_output_edge)->right);
                                }
                            tree_.left = position;
                            tree_.right = right;
                            // Must set return value before
                            // updating right, else the
                            // last tree will be skipped.
                            bool rv = current_input_edge < input_edge_end
                                      || position < maxpos;
                            position = right;
                            advanced_ = rv;
                        }
                    advanced_ = false;
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
