#ifndef FWDPP_TS_RECORDING_EDGE_BUFFER_HPP
#define FWDPP_TS_RECORDING_EDGE_BUFFER_HPP

#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <limits>
#include <cstdint>
#include <vector>
#include "../definitions.hpp"
#include <fwdpp/types/nested_forward_lists.hpp>

namespace fwdpp
{
    namespace ts
    {
        //constexpr std::int64_t EDGE_BUFFER_NULL = -1;

        struct birth_data
        {
            double left, right;
            table_index_t child;
            //std::int64_t next;

            birth_data(double l, double r, table_index_t c)
                : left{l}, right{r}, child{c}//, next{EDGE_BUFFER_NULL}
            {
            }
        };

        using edge_buffer = nested_forward_lists<birth_data, std::int32_t, -1>;

        // Below are functions for liftover of an edge buffer
        // to a table collection

        struct parent_location
        {
            table_index_t parent;
            std::size_t start, stop;
            parent_location(table_index_t p, std::size_t start_, std::size_t stop_)
                : parent{p}, start{start_}, stop{stop_}
            {
            }
        };

        template <typename TableCollectionType>
        inline std::vector<parent_location>
        find_pre_existing_edges(
            const TableCollectionType& tables,
            const std::vector<table_index_t>& alive_at_last_simplification,
            const fwdpp::ts::edge_buffer& new_edges)
        // FIXME: the indexing step need go no farther than the time of the most
        // recent node in alive_at_last_simplification.
        {
            std::vector<table_index_t> alive_with_new_edges;
            for (auto a : alive_at_last_simplification)
                {
                    if (new_edges.head(a) != edge_buffer::null)
                        {
                            alive_with_new_edges.push_back(a);
                        }
                }
            if (alive_with_new_edges.empty()) // get out early
                {
                    return {};
                }

            // index where each node already has edges.
            std::vector<std::size_t> starts(tables.num_nodes(),
                                            std::numeric_limits<std::size_t>::max()),
                stops(tables.num_nodes(), std::numeric_limits<std::size_t>::max());
            for (std::size_t i = 0; i < tables.num_edges(); ++i)
                {
                    if (starts[tables.edges[i].parent]
                        == std::numeric_limits<std::size_t>::max())
                        {
                            starts[tables.edges[i].parent] = i;
                            // FIXME: idiomatically, this should be i+1
                            stops[tables.edges[i].parent] = i;
                        }
                    else
                        { // FIXME: idiomatically, this should be i+1
                            stops[tables.edges[i].parent] = i;
                        }
                }

            std::vector<parent_location> existing_edges;
            for (auto a : alive_with_new_edges)
                {
                    existing_edges.emplace_back(a, starts[a], stops[a]);
                }

            // Our only sort!!
            std::sort(begin(existing_edges), end(existing_edges),
                      [&tables](const parent_location& lhs, const parent_location& rhs) {
                          // NOTE: have to take -time so that tuple sorting works
                          auto t0 = -tables.nodes[lhs.parent].time;
                          auto t1 = -tables.nodes[rhs.parent].time;
                          return std::tie(t0, lhs.start, lhs.parent)
                                 < std::tie(t1, rhs.start, rhs.parent);
                      });

            // FIXME: this should be debug only
            for (std::size_t i = 1; i < existing_edges.size(); ++i)
                {
                    auto t0 = tables.nodes[existing_edges[i - 1].parent].time;
                    auto t1 = tables.nodes[existing_edges[i].parent].time;
                    if (t0 < t1)
                        {
                            throw std::runtime_error(
                                "existing edges not properly sorted by time");
                        }
                }

            return existing_edges;
        }

        template <typename TableCollectionType>
        std::size_t
        handle_pre_existing_edges(
            const TableCollectionType& tables, const edge_buffer& new_edges,
            const std::vector<parent_location>& existing_edges,
            typename TableCollectionType::edge_table& edge_liftover)
        {
            std::size_t offset = 0;
            for (const auto& ex : existing_edges)
                {
                    // FIXME: this while loop is repeated 2x just w/different
                    // ranges
                    while (offset < tables.num_edges()
                           && tables.nodes[tables.edges[offset].parent].time
                                  > tables.nodes[ex.parent].time)
                        {
                            edge_liftover.emplace_back(tables.edges[offset]);
                            ++offset;
                        }
                    if (ex.start != std::numeric_limits<std::size_t>::max())
                        {
                            while (offset < ex.start
                                   && tables.nodes[tables.edges[offset].parent].time
                                          >= tables.nodes[ex.parent].time)
                                {
                                    edge_liftover.emplace_back(tables.edges[offset]);
                                    ++offset;
                                }
                            // FIXME: stop condition isn't idiomatic
                            for (decltype(ex.start) i = ex.start; i < ex.stop + 1; ++i)
                                {
                                    edge_liftover.emplace_back(tables.edges[i]);
                                }
                            offset = ex.stop + 1;
                        }
                    auto n = new_edges.head(ex.parent);
                    while (n != edge_buffer::null)
                        {
                            const auto & birth = new_edges.fetch(n);
                            edge_liftover.emplace_back(
                                typename TableCollectionType::edge_t{
                                    birth.left, birth.right,
                                    ex.parent, birth.child});
                            n = new_edges.next(n);
                        }
                }
            return offset;
        }

        template <typename TableCollectionType>
        inline void
        copy_births_since_last_simplification(
            const edge_buffer& new_edges, const TableCollectionType& tables,
            double max_time, typename TableCollectionType::edge_table& edge_liftover)
        {

            // TODO: should validate data in new_edges
            edge_liftover.clear(); // Re-use this buffer b/c it'll get big.

            // Go backwards through new births, and add them
            // to our temporary edge table if they are newer
            // than the last simplification time

            for (auto b = new_edges.rbegin(); b < new_edges.rend(); ++b)
                {
                    auto parent = new_edges.convert_to_head_index(b);
                    if(parent < 0)
                    {
                        throw std::runtime_error("negative parent value");
                    }
                    if(parent >= std::numeric_limits<table_index_t>::max())
                    {
                        throw std::overflow_error("parent value overflows");
                    }
                    auto cast_parent = static_cast<table_index_t>(parent);
                    auto ptime = tables.nodes[parent].time;
                    if (*b != edge_buffer::null && ptime > max_time)
                        {
                            auto n = *b;
                            while (n != edge_buffer::null)
                                {
                                    const auto & birth = new_edges.fetch(n);
                                    edge_liftover.emplace_back(
                                        typename TableCollectionType::edge_t{
                                            birth.left,
                                            birth.right, cast_parent,
                                            birth.child});
                                    n = new_edges.next(n);
                                }
                        }
                    else if (*b != edge_buffer::null && ptime <= max_time)
                        {
                            break;
                        }
                }
        }

        template <typename TableCollectionType>
        void
        stitch_together_edges(
            const std::vector<table_index_t>& alive_at_last_simplification,
            double max_time, edge_buffer& new_edges,
            typename TableCollectionType::edge_table& edge_liftover,
            TableCollectionType& tables)
        {
            copy_births_since_last_simplification(new_edges, tables, max_time,
                                                  edge_liftover);
            auto existing_edges = find_pre_existing_edges(
                tables, alive_at_last_simplification, new_edges);
            auto offset = handle_pre_existing_edges(tables, new_edges, existing_edges,
                                                    edge_liftover);
            for (; offset < tables.num_edges(); ++offset)
                {
                    edge_liftover.emplace_back(tables.edges[offset]);
                }
            tables.edges.assign(begin(edge_liftover), end(edge_liftover));
            // This resets sizes to 0, but keeps the memory allocated.
            edge_liftover.clear();
            new_edges.reset(tables.num_nodes());
        }
    }
}

#endif
