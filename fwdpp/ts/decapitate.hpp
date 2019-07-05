#ifndef FWDPP_TS_DECAPITATE_HPP
#define FWDPP_TS_DECAPITATE_HPP

#include <algorithm>
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ts
    {
        inline void
        decapitate(table_collection& tables, double time)
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
