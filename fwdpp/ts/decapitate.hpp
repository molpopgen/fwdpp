#ifndef FWDPP_TS_DECAPITATE_HPP
#define FWDPP_TS_DECAPITATE_HPP

#include <iostream>
#include <cmath>
#include <algorithm>
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        inline void
        decapitate(table_collection& tables, double time,
                   bool remove_mutations)
        /// \ brief "Decaptitate" a table_collection.
        ///
        /// Removes all nodes with time <= \a time
        /// as well as corresponding entries from
        /// the edge table.
        ///
        /// The table indexes are rebuilt.
        ///
        /// This function has little use outside of
        /// testing.
        {
            if (remove_mutations == true)
                {
                    auto itr = std::remove_if(
                        begin(tables.mutation_table),
                        end(tables.mutation_table),
                        [&tables, time](const mutation_record& mr) {
                            return tables.node_table[mr.node].time <= time;
                        });
                    tables.mutation_table.erase(itr,
                                                end(tables.mutation_table));
                    std::vector<int> site_referred_to(tables.site_table.size(),
                                                      0);
                    for (auto& mr : tables.mutation_table)
                        {
                            site_referred_to[mr.site]++;
                        }
                    for (std::size_t i = 0; i < site_referred_to.size(); ++i)
                        {
                            if (site_referred_to[i] == 0)
                                {
                                    tables.site_table[i].position
                                        = std::numeric_limits<
                                            double>::quiet_NaN();
                                }
                        }
                    auto site_itr = std::remove_if(
                        begin(tables.site_table), end(tables.site_table),
                        [](const site& s) { return std::isnan(s.position); });
                    tables.site_table.erase(site_itr, end(tables.site_table));
                }
            // NOTE: std::remove_if is stable w.r.to elements not removed
            tables.node_table.erase(std::remove_if(begin(tables.node_table),
                                                   end(tables.node_table),
                                                   [time](const node& n) {
                                                       return n.time <= time;
                                                   }),
                                    end(tables.node_table));
            tables.edge_table.erase(
                std::remove_if(begin(tables.edge_table),
                               end(tables.edge_table),
                               [&tables](const edge& e) {
                                   auto c = static_cast<std::size_t>(e.child)
                                            >= tables.node_table.size();
                                   auto p = static_cast<std::size_t>(e.parent)
                                            >= tables.node_table.size();
                                   return c || p;
                               }),
                end(tables.edge_table));
            tables.build_indexes();
        }
    } // namespace ts
} // namespace fwdpp

#endif
