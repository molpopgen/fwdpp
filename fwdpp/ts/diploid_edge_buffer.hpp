#ifndef FWDPP_TS_DIPLOID_EDGE_BUFFER_HPP
#define FWDPP_TS_DIPLOID_EDGE_BUFFER_HPP

#include <limits>
#include <tuple>
#include <stdexcept>
#include "edge_buffer.hpp"
#include "detail/split_breakpoints.hpp"

namespace fwdpp
{
    namespace ts
    {
        using diploid_edge_buffer = edge_buffer<2>;

        inline TS_NODE_INT
        register_offspring_and_buffer_edges(
            const std::vector<double>& recombination_breakpoints,
            std::size_t parent_index,
            std::tuple<TS_NODE_INT, TS_NODE_INT> parent_nodes,
            std::int32_t offspring_population, double birth_time,
            table_collection& tables, diploid_edge_buffer& buffer)
        {
            auto offspring_node
                = tables.emplace_back_node(offspring_population, birth_time);
            if (offspring_node
                >= std::numeric_limits<fwdpp::ts::TS_NODE_INT>::max())
                {
                    throw std::invalid_argument("node index too large");
                }
            std::size_t first_node_index = 0;
            auto first_node = std::get<0>(parent_nodes);
            auto second_node = std::get<1>(parent_nodes);
            if (first_node > second_node)
                {
                    first_node_index = 1;
                }
            detail::split_breakpoints(
                recombination_breakpoints, parent_nodes, offspring_node,
                tables.genome_length(),
                [parent_index, first_node_index, first_node,
                 &buffer](edge&& e) {
                    std::size_t offset = (e.parent == first_node)
                                             ? first_node_index
                                             : !first_node_index;
                    buffer.current_epoch[parent_index][offset].emplace_back(
                        std::move(e));
                });
            return offspring_node;
        }
    } // namespace ts
} // namespace fwdpp
#endif

