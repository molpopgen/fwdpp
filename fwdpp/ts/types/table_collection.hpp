#ifndef FWDPP_TS_DETAIL_TABLE_COLLECTION_HPP
#define FWDPP_TS_DETAIL_TABLE_COLLECTION_HPP

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
#include "node.hpp"
#include "edge.hpp"
#include "site.hpp"
#include "mutation_record.hpp"
#include "generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace types
        {
            template <typename SignedInteger> struct table_collection
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
                static constexpr SignedInteger null = generate_null_id<SignedInteger>();

                using id_type = SignedInteger;
                using node = types::node<SignedInteger>;
                using edge = types::edge<SignedInteger>;
                using site = types::site;
                using mutation_record = types::mutation_record<SignedInteger>;

                using node_table = std::vector<node>;
                using edge_table = std::vector<edge>;
                using site_table = std::vector<site>;
                using mutation_table = std::vector<mutation_record>;
                using edge_t [[deprecated("use ::edge instead")]] = edge;
                using node_t [[deprecated("use ::node instead")]] = node;
                using site_t [[deprecated("use ::site instead")]] = site;
                using mutation_t [[deprecated("use ::mutation_record instead")]]
                = mutation_record;

                /// Node table for this simulation
                node_table nodes;
                /// Edge table for this simulation
                edge_table edges;
                /// Mutation table for this simulation;
                mutation_table mutations;
                /// Site table
                site_table sites;
                /// The input edge vector. "I" in \cite Kelleher2016-cb, page 13
                std::vector<SignedInteger> input_left;
                /// The output edge vector. "O" in \cite Kelleher2016-cb, page 13
                std::vector<SignedInteger> output_right;
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

                table_collection(const id_type num_initial_nodes,
                                 const double initial_time, id_type pop,
                                 const double maxpos)
                    : L{maxpos}, nodes{}, edges{}, mutations{}, sites{}, input_left{},
                      output_right{}, edge_offset{0}
                {
                    if (maxpos <= 0 || !std::isfinite(maxpos))
                        {
                            throw std::invalid_argument("maxpos must be > 0 and finite");
                        }
                    for (id_type i = 0; i < num_initial_nodes; ++i)
                        {
                            nodes.push_back(node{pop, initial_time});
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

                id_type
                push_back_node(double time, std::int32_t pop)
                {
                    nodes.push_back(node{pop, time});
                    return static_cast<id_type>(nodes.size() - 1);
                }

                template <typename... args>
                id_type
                emplace_back_node(args&&... Args)
                {
                    nodes.emplace_back(node{std::forward<args>(Args)...});
                    return static_cast<id_type>(nodes.size() - 1);
                }

                std::size_t
                push_back_edge(double l, double r, id_type parent, id_type child)
                {
                    edges.push_back(edge{l, r, parent, child});
                    return edges.size();
                }

                template <typename... args>
                std::size_t
                emplace_back_edge(args&&... Args)
                {
                    edges.emplace_back(edge{std::forward<args>(Args)...});
                    return edges.size();
                }

                template <typename... args>
                std::size_t
                push_back_site(args&&... Args)
                {
                    sites.push_back(site{std::forward<args>(Args)...});
                    return sites.size() - 1;
                }

                template <typename... args>
                std::size_t
                emplace_back_site(args&&... Args)
                {
                    sites.emplace_back(site{std::forward<args>(Args)...});
                    return sites.size() - 1;
                }

                template <typename... args>
                std::size_t
                push_back_mutation(args&&... Args)
                {
                    mutations.push_back(mutation_record{std::forward<args>(Args)...});
                    return mutations.size() - 1;
                }

                template <typename... args>
                std::size_t
                emplace_back_mutation(args&&... Args)
                {
                    mutations.emplace_back(mutation_record{std::forward<args>(Args)...});
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
                    std::sort(begin(input_left), end(input_left),
                              [this](auto i, auto j) {
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
        } // namespace detail
    }     // namespace ts
} // namespace fwdpp

#endif
