#ifndef FWDPP_TS_DECAPITATE_HPP
#define FWDPP_TS_DECAPITATE_HPP

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <vector>

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType>
        inline void
        decapitate(TableCollectionType& tables, double time, bool remove_mutations)
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
                        begin(tables.mutations), end(tables.mutations),
                        [&tables, time](const typename TableCollectionType::
                                            mutation_table::value_type& mr) {
                            return tables.nodes[mr.node].time <= time;
                        });
                    tables.mutations.erase(itr, end(tables.mutations));
                    std::vector<int> site_referred_to(tables.sites.size(), 0);
                    for (auto& mr : tables.mutations)
                        {
                            site_referred_to[mr.site]++;
                        }
                    for (std::size_t i = 0; i < site_referred_to.size(); ++i)
                        {
                            if (site_referred_to[i] == 0)
                                {
                                    tables.sites[i].position
                                        = std::numeric_limits<double>::quiet_NaN();
                                }
                        }
                    auto site_itr = std::remove_if(
                        begin(tables.sites), end(tables.sites),
                        [](const typename TableCollectionType::site_table::value_type&
                               s) { return std::isnan(s.position); });
                    tables.sites.erase(site_itr, end(tables.sites));
                }
            // NOTE: std::remove_if is stable w.r.to elements not removed
            tables.nodes.erase(
                std::remove_if(
                    begin(tables.nodes), end(tables.nodes),
                    [time](
                        const typename TableCollectionType::node_table::value_type& n) {
                        return n.time <= time;
                    }),
                end(tables.nodes));
            tables.edges.erase(
                std::remove_if(
                    begin(tables.edges), end(tables.edges),
                    [&tables](
                        const typename TableCollectionType::edge_table::value_type& e) {
                        auto c
                            = static_cast<std::size_t>(e.child) >= tables.nodes.size();
                        auto p
                            = static_cast<std::size_t>(e.parent) >= tables.nodes.size();
                        return c || p;
                    }),
                end(tables.edges));
            tables.build_indexes();
        }
    } // namespace ts
} // namespace fwdpp

#endif
