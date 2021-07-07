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
            return [&tables](const typename TableCollectionType::edge& a,
                             const typename TableCollectionType::edge& b) {
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
            return [&tables](const typename TableCollectionType::edge& a,
                             const typename TableCollectionType::edge& b) {
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
        record_site_during_rebuild(const typename TableCollectionType::site& s,
                                   typename TableCollectionType::mutation_record& mr,
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

        template <typename TableCollectionType>
        std::size_t
        count_trees(const TableCollectionType& tables)
        /// Count the number of trees.
        /// Tables must be indexed.
        ///
        /// Implementation ported from forrustts.
        {
            if (!tables.indexed())
                {
                    throw tables_error("tables are not indexed");
                }
            std::size_t num_trees = 0;

            auto input_edge = begin(tables.input_left);
            auto output_edge = begin(tables.output_right);
            const auto input_end = end(tables.input_left);
            const auto output_end = end(tables.output_right);

            double tree_left = 0.0;

            while (input_edge < input_end || tree_left < tables.genome_length())
                {
                    while (output_edge < output_end
                           && tables.edges[*output_edge].right == tree_left)
                        {
                            ++output_edge;
                        }
                    while (input_edge < input_end
                           && tables.edges[*input_edge].left == tree_left)
                        {
                            ++input_edge;
                        }
                    auto tree_right = tables.genome_length();
                    if (input_edge < input_end)
                        {
                            tree_right
                                = std::min(tree_right, tables.edges[*input_edge].left);
                        }
                    if (output_edge < output_end)
                        {
                            tree_right
                                = std::min(tree_right, tables.edges[*output_edge].right);
                        }
                    tree_left = tree_right;
                    ++num_trees;
                }
            return num_trees;
        }
    }
}

#endif
