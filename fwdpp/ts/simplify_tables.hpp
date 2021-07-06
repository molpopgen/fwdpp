#ifndef FWDPP_TS_SIMPLIFY_TABLES_HPP__
#define FWDPP_TS_SIMPLIFY_TABLES_HPP__

#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "definitions.hpp"
#include "recording/edge_buffer.hpp"
#include "simplification/simplification.hpp"
#include "types/table_collection.hpp"
#include <fwdpp/ts/simplification_flags.hpp>
#include <fwdpp/ts/simplify_tables_output.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename SignedInteger>
        inline void
        simplify_tables(const std::vector<SignedInteger>& samples,
                        const simplification_flags /*flags*/,
                        simplification::simplifier_internal_state<SignedInteger>& state,
                        types::table_collection<SignedInteger>& input_tables,
                        simplify_tables_output<SignedInteger>& simplification_output)
        {
            static_assert(std::is_integral<SignedInteger>::value,
                          "SignedInteger must be an integer type");
            static_assert(std::is_signed<SignedInteger>::value,
                          "SignedInteger must be a signed type");
            state.clear();
            state.ancestry.reset(input_tables.nodes.size());
            simplification_output.clear();
            simplification_output.idmap.resize(input_tables.nodes.size());
            std::fill(begin(simplification_output.idmap),
                      end(simplification_output.idmap),
                      types::table_collection<SignedInteger>::null);

            // We take our samples and add them to both the output
            // node list and initialize their ancestry with
            // a segment on [0,L).
            simplification::record_sample_nodes(samples, input_tables, state,
                                                simplification_output.idmap);
            // At this point, our edges are sorted by birth
            // order of parents, from present to past.
            // We can now work our way up the pedigree.
            // This outer loop differs from how we describe it in the
            // paper, but the strict sorting of edges means that this
            // equivalent.
            auto edge_ptr = input_tables.edges.cbegin();
            const auto edge_end = input_tables.edges.cend();
            auto new_edge_destination = begin(input_tables.edges);
            while (edge_ptr < edge_end)
                {
                    auto u = edge_ptr->parent;
                    edge_ptr = simplification::find_parent_child_segment_overlap(
                        input_tables.genome_length(), edge_ptr, edge_end, u, state);
                    simplification::merge_ancestors(input_tables.genome_length(),
                                                    input_tables.nodes, u, state,
                                                    simplification_output.idmap);
                    if (state.new_edge_table.size() >= 1024
                        && new_edge_destination + state.new_edge_table.size() < edge_ptr)
                        {
                            new_edge_destination = std::copy(begin(state.new_edge_table),
                                                             end(state.new_edge_table),
                                                             new_edge_destination);
                            state.new_edge_table.clear();
                        }
                }
            if (state.new_edge_table.empty() == false)
                {
                    new_edge_destination
                        = std::copy(begin(state.new_edge_table),
                                    end(state.new_edge_table), new_edge_destination);
                }

            // FIXME: figure out what to do with this?
            assert(
                static_cast<std::size_t>(std::count_if(
                    begin(simplification_output.idmap), end(simplification_output.idmap),
                    [](const SignedInteger i) {
                        return i != types::table_collection<SignedInteger>::null;
                    }))
                == state.new_node_table.size());

            simplification::simplify_mutations(
                state, input_tables, simplification_output.preserved_mutations);
            simplification::transfer_new_nodes_and_edges(new_edge_destination, state,
                                                         input_tables);
        }

        template <typename SignedInteger>
        inline void
        simplify_tables(const std::vector<SignedInteger>& samples,
                        const std::vector<SignedInteger>& alive_at_last_simplification,
                        simplification_flags /*flags*/,
                        simplification::simplifier_internal_state<SignedInteger>& state,
                        types::table_collection<SignedInteger>& input_tables,
                        edge_buffer<SignedInteger>& buffer,
                        simplify_tables_output<SignedInteger>& simplification_output)
        {
            state.clear();
            state.ancestry.reset(input_tables.nodes.size());
            simplification_output.clear();
            simplification_output.idmap.resize(input_tables.nodes.size());
            std::fill(begin(simplification_output.idmap),
                      end(simplification_output.idmap),
                      types::table_collection<SignedInteger>::null);

            // We take our samples and add them to both the output
            // node list and initialize their ancestry with
            // a segment on [0,L).
            simplification::record_sample_nodes(samples, input_tables, state,
                                                simplification_output.idmap);

            // NOTE: it may be better to require that alive_at_last_simplification
            // not be empty, so that we force it to be filled with
            // the initial nodes for the first simplification, too.

            // Process all edges since the last simplification.
            double max_time = -std::numeric_limits<double>::max();
            for (auto a : alive_at_last_simplification)
                {
                    max_time = std::max(max_time, input_tables.nodes[a].time);
                }
            auto buffer_rend = buffer.rbegin();
            for (; buffer_rend < buffer.rend(); ++buffer_rend)
                {
                    auto parent = buffer.convert_to_head_index(buffer_rend);
                    auto ptime = input_tables.nodes[parent].time;
                    if (*buffer_rend != edge_buffer<SignedInteger>::null
                        && ptime > max_time)
                        // Then *b is a parent node born after the last
                        // simplification that did leave offspring
                        {
                            state.overlapper.clear_queue();
                            auto n = *buffer_rend;
                            simplification::process_births_from_buffer(
                                n, buffer, state.ancestry, state.overlapper);
                            state.overlapper.finalize_queue(
                                input_tables.genome_length());

                            simplification::merge_ancestors(
                                input_tables.genome_length(), input_tables.nodes,
                                static_cast<SignedInteger>(parent), state,
                                simplification_output.idmap);
                        }
                    else if (*buffer_rend != edge_buffer<SignedInteger>::null
                             && ptime <= max_time)
                        {
                            break;
                        }
                }

            // Now comes the tricky bit, which is to handle edges whose
            // parent nodes were alive at the time of the last simplification.

            auto existing_edges = find_pre_existing_edges(
                input_tables, alive_at_last_simplification, buffer);
            auto edge_ptr = input_tables.edges.cbegin();
            const auto edge_end = input_tables.edges.cend();
            for (auto& ex : existing_edges)
                {
                    while (edge_ptr < edge_end
                           && input_tables.nodes[edge_ptr->parent].time
                                  > input_tables.nodes[ex.parent].time)
                        {
                            auto u = edge_ptr->parent;
                            edge_ptr = simplification::find_parent_child_segment_overlap(
                                input_tables.genome_length(), edge_ptr, edge_end, u,
                                state);
                            simplification::merge_ancestors(input_tables.genome_length(),
                                                            input_tables.nodes, u, state,
                                                            simplification_output.idmap);
                        }
                    if (ex.start != std::numeric_limits<std::size_t>::max())
                        {
                            // FIXME: maybe we should just use ptrdiff_t for ex.start/stop?
                            auto offset = static_cast<std::ptrdiff_t>(ex.start);
                            auto end_of_range = input_tables.edges.cbegin() + offset;
                            while (edge_ptr < end_of_range
                                   // FIXME: the next parts of this check should
                                   // not be necessary
                                   && input_tables.nodes[edge_ptr->parent].time
                                          >= input_tables.nodes[ex.parent].time)
                                {
                                    auto u = edge_ptr->parent;
                                    edge_ptr = simplification::
                                        find_parent_child_segment_overlap(
                                            input_tables.genome_length(), edge_ptr,
                                            edge_end, u, state);
                                    simplification::merge_ancestors(
                                        input_tables.genome_length(), input_tables.nodes,
                                        u, state, simplification_output.idmap);
                                }
                            if (edge_ptr != input_tables.edges.cbegin() + offset)
                                {
                                    throw std::runtime_error(
                                        "unexpected location in input edge table");
                                }
                            if (edge_ptr->parent != ex.parent)
                                {
                                    throw std::runtime_error(
                                        "unexpected parent value in input edge table");
                                }
                        }
                    // now handle ex.parent
                    state.overlapper.clear_queue();
                    if (ex.start != std::numeric_limits<std::size_t>::max())
                        {
                            auto offset = static_cast<std::ptrdiff_t>(ex.stop + 1);
                            auto end_of_range = input_tables.edges.cbegin() + offset;
                            for (; edge_ptr < end_of_range; ++edge_ptr)
                                {
                                    if (edge_ptr->parent != ex.parent)
                                        {
                                            throw std::runtime_error(
                                                "processing unexpected parent node");
                                        }
                                    simplification::queue_children(
                                        edge_ptr->child, edge_ptr->left, edge_ptr->right,
                                        state.ancestry, state.overlapper);
                                }
                            if (edge_ptr < edge_end && edge_ptr->parent == ex.parent)
                                {
                                    throw std::runtime_error(
                                        "error traversing pre-existing edges for "
                                        "parent");
                                }
                        }
                    auto n = buffer.head(ex.parent);
                    simplification::process_births_from_buffer(n, buffer, state.ancestry,
                                                               state.overlapper);
                    state.overlapper.finalize_queue(input_tables.genome_length());
                    simplification::merge_ancestors(input_tables.genome_length(),
                                                    input_tables.nodes, ex.parent, state,
                                                    simplification_output.idmap);
                }

            // Handle the remaning edges
            while (edge_ptr < edge_end)
                {
                    auto u = edge_ptr->parent;
                    edge_ptr = simplification::find_parent_child_segment_overlap(
                        input_tables.genome_length(), edge_ptr, edge_end, u, state);
                    simplification::merge_ancestors(input_tables.genome_length(),
                                                    input_tables.nodes, u, state,
                                                    simplification_output.idmap);
                }

            // FIXME: figure out what to do with this?
            assert(
                static_cast<std::size_t>(std::count_if(
                    begin(simplification_output.idmap), end(simplification_output.idmap),
                    [](const SignedInteger i) {
                        return i != types::table_collection<SignedInteger>::null;
                    }))
                == state.new_node_table.size());

            simplification::simplify_mutations(
                state, input_tables, simplification_output.preserved_mutations);
            input_tables.edges.resize(state.new_edge_table.size());
            std::move(begin(state.new_edge_table), end(state.new_edge_table),
                      begin(input_tables.edges));
            input_tables.nodes.resize(state.new_node_table.size());
            std::move(begin(state.new_node_table), end(state.new_node_table),
                      begin(input_tables.nodes));
            buffer.reset(input_tables.num_nodes());
        }

        // Below, we have deprecated functions that take multiple outputs
        // as arguments.  In 0.10.0, we introduced simplify_tables_output_t
        // to hold all of these in a single object.

        template <typename SignedInteger, typename PreservedVariantIndexes>
        [[deprecated]] inline void
        simplify_tables(const SignedInteger& samples, const simplification_flags flags,
                        simplification::simplifier_internal_state<SignedInteger>& state,
                        types::table_collection<SignedInteger>& input_tables,
                        std::vector<SignedInteger>& idmap,
                        PreservedVariantIndexes& preserved_variants)
        {
            simplify_tables_output_t<SignedInteger, PreservedVariantIndexes>
                simplification_output;
            simplify_tables(samples, flags, state, input_tables, simplification_output);
            idmap.swap(simplification_output.idmap);
            preserved_variants.swap(simplification_output.preserved_mutations);
        }

        template <typename SignedInteger, typename PreservedVariantIndexes>
        [[deprecated]] inline void
        simplify_tables(const std::vector<SignedInteger>& samples,
                        types::table_collection<SignedInteger>& input_tables,
                        simplification_flags flags, std::vector<SignedInteger>& idmap,
                        PreservedVariantIndexes& preserved_variants)
        {
            auto state = simplification::make_simplifier_internal_state(input_tables);
            simplify_tables_output_t<SignedInteger, PreservedVariantIndexes>
                simplification_output;
            simplify_tables(samples, flags, state, input_tables, simplification_output);
            idmap.swap(simplification_output.idmap);
            preserved_variants.swap(simplification_output.preserved_mutations);
        }

        template <typename SignedInteger>
        [[deprecated]] inline void
        simplify_tables(const std::vector<SignedInteger>& samples,
                        const simplification_flags flags,
                        types::table_collection<SignedInteger>& input_tables)
        {
            std::vector<SignedInteger> idmap;
            simplify_tables_output_t<SignedInteger, std::vector<std::size_t>>
                simplification_output;
            auto state = simplification::make_simplifier_internal_state(input_tables);
            simplify_tables(samples, flags, state, input_tables, simplification_output);
            idmap.swap(simplification_output.idmap);
        }

        template <typename SignedInteger, typename PreservedVariantIndexes>
        [[deprecated]] inline void
        simplify_tables(const std::vector<SignedInteger>& samples,
                        const std::vector<SignedInteger>& alive_at_last_simplification,
                        simplification_flags flags,
                        simplification::simplifier_internal_state<SignedInteger>& state,
                        types::table_collection<SignedInteger>& input_tables,
                        edge_buffer<SignedInteger>& buffer,
                        std::vector<SignedInteger>& idmap,
                        PreservedVariantIndexes& preserved_variants)
        {
            simplify_tables_output_t<SignedInteger, PreservedVariantIndexes>
                simplification_output;
            simplify_tables(samples, alive_at_last_simplification, flags, state,
                            input_tables, buffer, simplification_output);
            idmap.swap(simplification_output.idmap);
            preserved_variants.swap(simplification_output.preserved_mutations);
        }

    }
}

#endif
