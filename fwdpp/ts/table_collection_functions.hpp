#ifndef FWDPP_TS_TABLE_COLLECTION_FUNCTIONS_HPP
#define FWDPP_TS_TABLE_COLLECTION_FUNCTIONS_HPP

#include <tuple>
#include <cstddef>
#include <algorithm>
#include <stdexcept>
#include <fwdpp/ts/exceptions.hpp>

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType>
        inline auto
        get_edge_sort_cmp(TableCollectionType& tables) noexcept
        /*!
         * @version 0.9.0 Added to library
         */
        {
            return [&tables](const typename TableCollectionType::edge_t& a,
                             const typename TableCollectionType::edge_t& b) {
                auto ga = tables.nodes[a.parent].time;
                auto gb = tables.nodes[b.parent].time;
                if (ga == gb)
                    {
                        if (a.parent == b.parent)
                            {
                                if (a.child == b.child)
                                    {
                                        return a.left < b.left;
                                    }
                                return a.child < b.child;
                            }
                        return a.parent < b.parent;
                    }
                return ga > gb;
            };
        }

        template <typename TableCollectionType>
        inline auto
        get_minimal_edge_sort_cmp(const TableCollectionType& tables)
        {
            return [&tables](const typename TableCollectionType::edge_t& a,
                             const typename TableCollectionType::edge_t& b) {
                auto ga = tables.nodes[a.parent].time;
                auto gb = tables.nodes[b.parent].time;
                return ga > gb
                       && (std::tie(a.child, a.left) < std::tie(b.child, b.left));
            };
        }

        template <typename TableCollectionType>
        inline void
        sort_edge_table(std::ptrdiff_t offset, TableCollectionType& tables)
        {
            if (offset < 0 || offset >= static_cast<std::ptrdiff_t>(tables.edges.size()))
                {
                    throw std::out_of_range("invalid edge table offset");
                }
            auto cmp = get_edge_sort_cmp(tables);
            std::sort(tables.edges.begin() + offset, tables.edges.end(), cmp);
            if (offset > 0)
                {
                    std::rotate(begin(tables.edges), begin(tables.edges) + offset,
                                end(tables.edges));
                }
        }

        template <typename TableCollectionType>
        inline void
        sort_edge_table(TableCollectionType& tables) noexcept
        {
            sort_edge_table(0, tables);
        }

        template <typename TableCollectionType>
        inline bool
        edge_table_strictly_sorted(const TableCollectionType& tables)
        {
            auto cmp = get_edge_sort_cmp(tables);
            return std::is_sorted(begin(tables.edges), end(tables.edges), cmp);
        }

        template <typename TableCollectionType>
        inline bool
        edge_table_minimally_sorted(const TableCollectionType& tables)
        {
            auto cmp = get_minimal_edge_sort_cmp(tables);
            return std::is_sorted(begin(tables.edges), end(tables.edges), cmp);
        }

        template <typename TableCollectionType>
        inline void
        sort_mutation_table(TableCollectionType& tables)
        {
            //mutations are sorted by increasing position
            std::sort(begin(tables.mutations), end(tables.mutations),
                      [&tables](const auto& a, const auto& b) {
                          return tables.sites[a.site].position
                                 < tables.sites[b.site].position;
                      });
        }

        template <typename TableCollectionType>
        inline void
        record_site_during_rebuild(const typename TableCollectionType::site_t& s,
                                   typename TableCollectionType::mutation_t& mr,
                                   TableCollectionType& tables)
        {
            if (tables.sites.empty() || tables.sites.back().position != s.position)
                {
                    tables.sites.push_back(s);
                }
            mr.site = tables.sites.size() - 1;
        }

        template <typename TableCollectionType>
        inline void
        rebuild_site_table(TableCollectionType& tables)
        /// Complete rebuild of the site table.
        {
            auto site_table_copy(tables.sites);
            tables.sites.clear();
            for (auto& mr : tables.mutations)
                {
                    auto os = mr.site;
                    record_site_during_rebuild(site_table_copy[mr.site], mr, tables);
                    if (site_table_copy[os].position != tables.sites[mr.site].position)
                        {
                            throw tables_error("error rebuilding site table");
                        }
                }
        }

        template <typename TableCollectionType>
        inline void
        sort_mutation_table_and_rebuild_site_table(TableCollectionType& tables)
        /// O(n*log(n)) time plus O(n) additional memory.
        /// This is called by sort_tables, but also should
        /// be called after adding mutations by some manual
        /// means to a table collection.
        {
            sort_mutation_table(tables);
            rebuild_site_table(tables);
        }

        template <typename TableCollectionType>
        inline void
        sort_tables_for_simplification(std::ptrdiff_t edge_table_offset,
                                       TableCollectionType& tables)
        /// Sorts the tables for simplification, which means only
        /// sorting edge and mutation tables, as the site table
        /// will be rebuilt during simplification.
        ///
        /// If \a edge_table_offset < 0, the edge table is not sorted.
        {
            if (edge_table_offset >= 0)
                {
                    sort_edge_table(edge_table_offset, tables);
                }
            sort_mutation_table(tables);
        }

        template <typename TableCollectionType>
        inline void
        sort_tables(std::ptrdiff_t edge_table_offset, TableCollectionType& tables)
        /// Sort all tables.  The site table is rebuilt.
        {
            sort_edge_table(edge_table_offset, tables);
            sort_mutation_table_and_rebuild_site_table(tables);
        }
    }
}

#endif
