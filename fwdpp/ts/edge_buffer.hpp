#ifndef FWDPP_TS_EDGE_BUFFER_HPP
#define FWDPP_TS_EDGE_BUFFER_HPP

#include <vector>
#include <utility>
#include <array>
#include "table_types.hpp"
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <std::size_t PLOIDY> struct edge_buffer
        {
            std::vector<std::array<edge_vector, PLOIDY>> current_epoch;
            edge_vector buffer;
            std::vector<std::pair<std::size_t, std::size_t>> offsets;

            edge_buffer() : current_epoch{}, buffer{}, offsets{} {}

            inline void
            begin_epoch(std::size_t n)
            {
                current_epoch.resize(n);
            }

            inline void
            end_epoch()
            {
                std::size_t epoch_offset = buffer.size();
                for (auto&& i : current_epoch)
                    {
                        for (auto&& j : i)
                            {
                                buffer.insert(end(buffer), begin(j), end(j));
                                j.clear();
                            }
                    }
                offsets.emplace_back(epoch_offset, buffer.size());
            };

            inline void
            prepare_tables_for_simplification(table_collection& tables)
            {
                for (auto i = offsets.rbegin(); i < offsets.rend(); ++i)
                    {
                        tables.edge_table.insert(end(tables.edge_table),
                                                 begin(buffer) + i->first,
                                                 begin(buffer) + i->second);
                    }
                buffer.clear();
                offsets.clear();

                tables.sort_mutations();
                if (tables.edge_offset > 0)
                    {
                        std::rotate(begin(tables.edge_table),
                                    begin(tables.edge_table)
                                        + tables.edge_offset,
                                    end(tables.edge_table));
                    }
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif

