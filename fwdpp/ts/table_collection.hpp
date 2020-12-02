#ifndef FWDPP_TS_TABLE_COLLECTION_HPP
#define FWDPP_TS_TABLE_COLLECTION_HPP

#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/ts/exceptions.hpp>
#include "definitions.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename NodeTableType, typename EdgeTableType, typename SiteTableType,
                  typename MutationTableType>
        struct table_collection
        /*!
		 * \brief A collection of tables for a single simulation.
         *
         * \version 0.7.0 Added to fwdpp
         * \version 0.9.0 Made a template class
		 */
        {
          private:
            /// Length of the genomic region.
            double L;

          public:
            using node_table = NodeTableType;
            using edge_table = EdgeTableType;
            using site_table = SiteTableType;
            using mutation_table = MutationTableType;
            using edge_t = typename EdgeTableType::value_type;
            using node_t = typename NodeTableType::value_type;
            using site_t = typename SiteTableType::value_type;
            using mutation_t = typename MutationTableType::value_type;

            /// Node table for this simulation
            node_table nodes;
            /// Edge table for this simulation
            edge_table edges;
            /// Mutation table for this simulation;
            mutation_table mutations;
            /// Site table
            site_table sites;
            /// The input edge vector. "I" in \cite Kelleher2016-cb, page 13
            std::vector<table_index_t> input_left;
            /// The output edge vector. "O" in \cite Kelleher2016-cb, page 13
            std::vector<table_index_t> output_right;
            /// This reflects the length of
            /// tables.edges after last simplification.
            /// It can be used to make sure we only sort
            /// newly-added nodes.
            std::ptrdiff_t edge_offset;

            explicit table_collection(const double maxpos)
                : L{maxpos}, nodes{}, edges{}, mutations{}, sites{}, input_left{},
                  output_right{}, edge_offset{0}
            {
                if (maxpos <= 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument("maxpos must be > 0 and finite");
                    }
            }

            table_collection(const table_index_t num_initial_nodes,
                             const double initial_time, table_index_t pop,
                             const double maxpos)
                : L{maxpos}, nodes{}, edges{}, mutations{}, sites{}, input_left{},
                  output_right{}, edge_offset{0}
            {
                if (maxpos <= 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument("maxpos must be > 0 and finite");
                    }
                for (table_index_t i = 0; i < num_initial_nodes; ++i)
                    {
                        nodes.push_back(node_t{pop, initial_time});
                    }
            }

            void
            clear() noexcept
            /// Clears internal vectors.
            {
                nodes.clear();
                edges.clear();
                mutations.clear();
                sites.clear();
            }

            table_index_t
            push_back_node(double time, std::int32_t pop)
            {
                nodes.push_back(node_t{pop, time});
                return static_cast<table_index_t>(nodes.size() - 1);
            }

            template <typename... args>
            table_index_t
            emplace_back_node(args&&... Args)
            {
                nodes.emplace_back(node_t{std::forward<args>(Args)...});
                return static_cast<table_index_t>(nodes.size() - 1);
            }

            std::size_t
            push_back_edge(double l, double r, table_index_t parent, table_index_t child)
            {
                edges.push_back(edge_t{l, r, parent, child});
                return edges.size();
            }

            template <typename... args>
            std::size_t
            emplace_back_edge(args&&... Args)
            {
                edges.emplace_back(edge_t{std::forward<args>(Args)...});
                return edges.size();
            }

            template <typename... args>
            std::size_t
            push_back_site(args&&... Args)
            {
                sites.push_back(site_t{std::forward<args>(Args)...});
                return sites.size() - 1;
            }

            template <typename... args>
            std::size_t
            emplace_back_site(args&&... Args)
            {
                sites.emplace_back(site_t{std::forward<args>(Args)...});
                return sites.size() - 1;
            }

            template <typename... args>
            std::size_t
            push_back_mutation(args&&... Args)
            {
                mutations.push_back(mutation_t{std::forward<args>(Args)...});
                return mutations.size() - 1;
            }

            template <typename... args>
            std::size_t
            emplace_back_mutation(args&&... Args)
            {
                mutations.emplace_back(mutation_t{std::forward<args>(Args)...});
                return mutations.size() - 1;
            }

            void
            build_indexes()
            /// Generates the index vectors referred to
            /// as I and O in Kelleher et al. (2016)
            {
                input_left.clear();
                output_right.clear();
                input_left.resize(edges.size());
                output_right.resize(edges.size());
                std::iota(begin(input_left), end(input_left), 0);
                std::iota(begin(output_right), end(output_right), 0);
                std::sort(begin(input_left), end(input_left), [this](auto i, auto j) {
                    if (edges[i].left == edges[j].left)
                        {
                            return nodes[edges[i].parent].time
                                   > nodes[edges[j].parent].time;
                        }
                    return edges[i].left < edges[j].left;
                });
                std::sort(begin(output_right), end(output_right),
                          [this](auto i, auto j) {
                              if (edges[i].right == edges[j].right)
                                  {
                                      return nodes[edges[i].parent].time
                                             < nodes[edges[j].parent].time;
                                  }
                              return edges[i].right < edges[j].right;
                          });
            }

            std::size_t
            num_nodes() const
            {
                return nodes.size();
            }

            std::size_t
            num_edges() const
            {
                return edges.size();
            }

            double
            genome_length() const
            {
                return L;
            }

            inline bool
            operator==(const table_collection& b) const
            {
                return genome_length() == b.genome_length() && sites == b.sites
                       && edges == b.edges && nodes == b.nodes
                       && mutations == b.mutations;
            }

            inline bool
            operator!=(const table_collection& b) const
            {
                return !(*this == b);
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
