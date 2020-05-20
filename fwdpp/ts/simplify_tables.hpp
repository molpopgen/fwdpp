#ifndef FWDPP_TS_SIMPLIFY_TABLES_HPP__
#define FWDPP_TS_SIMPLIFY_TABLES_HPP__

#include <cstdint>
#include <type_traits>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "definitions.hpp"
#include "simplification/simplification.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType, typename NodeVector,
                  typename PreservedVariantIndexes>
        inline void
        simplify_tables(
            const NodeVector& samples,
            simplification::simplifier_internal_state<TableCollectionType>& state,
            TableCollectionType& input_tables, NodeVector& idmap,
            PreservedVariantIndexes& preserved_variants)
        {
            static_assert(std::is_integral<typename NodeVector::value_type>::value,
                          "NodeVector::value type must be an integer type");
            static_assert(std::is_signed<typename NodeVector::value_type>::value,
                          "NodeVector::value type must be a signed type");
            state.clear();
            state.ancestry.init(input_tables.nodes.size());
            idmap.resize(input_tables.nodes.size());
            std::fill(begin(idmap), end(idmap), TS_NULL_NODE);

            // We take our samples and add them to both the output
            // node list and initialize their ancestry with
            // a segment on [0,L).
            simplification::record_sample_nodes(samples, input_tables, state, idmap);
            // Add samples for any preserved nodes in the tables:
            simplification::record_sample_nodes(input_tables.preserved_nodes,
                                                input_tables, state, idmap);
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
                                                    input_tables.nodes, u, state, idmap);
                    new_edge_destination = simplification::transfer_simplified_edges(
                        new_edge_destination, edge_ptr, 1024ul, state);
                }
            new_edge_destination = simplification::transfer_simplified_edges(
                new_edge_destination, edge_ptr, 0ul, state);

            // When there are preserved nodes, we need to re map
            // their input ids to output ids
            // FIXME/NOTE: should we be doing this here,
            // or expect the caller to do it?
            for (auto& p : input_tables.preserved_nodes)
                {
                    if (idmap[p] == TS_NULL_NODE)
                        {
                            throw std::runtime_error(
                                "preserved node output id maps to null");
                        }
                    p = idmap[p];
                }

            // FIXME: figure out what to do with this?
            assert(static_cast<std::size_t>(std::count_if(
                       begin(idmap), end(idmap),
                       [](const TS_NODE_INT i) { return i != TS_NULL_NODE; }))
                   == state.new_node_table.size());

            simplification::simplify_mutations(state, input_tables, preserved_variants);
            simplification::transfer_new_nodes_and_edges(new_edge_destination, state,
                                                         input_tables);
        }

        template <typename TableCollectionType, typename NodeVector,
                  typename PreservedVariantIndexes>
        inline void
        simplify_tables(const NodeVector& samples, TableCollectionType& input_tables,
                        NodeVector& idmap, PreservedVariantIndexes& preserved_variants)
        {
            auto state = simplification::make_simplifier_internal_state(input_tables);
            simplify_tables(samples, state, input_tables, idmap, preserved_variants);
        }

        template <typename TableCollectionType, typename NodeVector>
        inline void
        simplify_tables(const NodeVector& samples, TableCollectionType& input_tables)
        {
            NodeVector idmap;
            std::vector<std::size_t> preserved_variants;
            simplify_tables(samples, input_tables, idmap, preserved_variants);
        }

    }
}

#endif
