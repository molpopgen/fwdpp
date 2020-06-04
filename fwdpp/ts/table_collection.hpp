#ifndef FWDPP_TS_TABLE_COLLECTION_HPP
#define FWDPP_TS_TABLE_COLLECTION_HPP

#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/ts/exceptions.hpp>
#include "definitions.hpp"
#include "indexed_edge.hpp" //TODO: create fewer header dependencies

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
            void
            split_breakpoints_add_edges(
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const TS_NODE_INT next_index)
            {
                if (breakpoints.front() != 0.0)
                    {
                        this->push_back_edge(0., breakpoints.front(),
                                             std::get<0>(parents), next_index);
                    }
                // TODO: replace with exception via a debug mode
                assert(std::count(begin(breakpoints), end(breakpoints),
                                  std::numeric_limits<double>::max())
                       == 1);
                assert(breakpoints.back() == std::numeric_limits<double>::max());
                for (unsigned j = 1; j < breakpoints.size(); ++j)
                    {
                        double a = breakpoints[j - 1];
                        double b = (j < breakpoints.size() - 1) ? breakpoints[j] : L;
                        if (b <= a)
                            {
                                throw std::invalid_argument("right must be > left");
                            }
                        if (j % 2 == 0.)
                            {
                                this->push_back_edge(a, b, std::get<0>(parents),
                                                     next_index);
                            }
                        else
                            {
                                this->push_back_edge(a, b, std::get<1>(parents),
                                                     next_index);
                            }
                    }
            }

            void
            split_breakpoints(const std::vector<double>& breakpoints,
                              const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                              const TS_NODE_INT next_index)
            {
                if (breakpoints.empty())
                    {
                        this->push_back_edge(0., L, std::get<0>(parents), next_index);
                        return;
                    }
                auto itr
                    = std::adjacent_find(std::begin(breakpoints), std::end(breakpoints));
                if (itr == std::end(breakpoints))
                    {
                        split_breakpoints_add_edges(breakpoints, parents, next_index);
                    }
                else
                    {
                        // Here, we need to reduce the input
                        // breakpoints to only those seen
                        // an odd number of times.
                        // Even numbers of the same breakpoint
                        // are "double x-overs" and thus
                        // cannot affect the genealogy.
                        std::vector<double> odd_breakpoints;
                        auto start = breakpoints.begin();
                        while (itr < breakpoints.end())
                            {
                                auto not_equal = std::find_if(
                                    itr, breakpoints.end(),
                                    [itr](const double d) { return d != *itr; });
                                int even = (std::distance(itr, not_equal) % 2 == 0.0);
                                odd_breakpoints.insert(odd_breakpoints.end(), start,
                                                       itr + 1 - even);
                                start = not_equal;
                                itr = std::adjacent_find(start, std::end(breakpoints));
                            }
                        odd_breakpoints.insert(odd_breakpoints.end(), start,
                                               breakpoints.end());
                        split_breakpoints_add_edges(odd_breakpoints, parents,
                                                    next_index);
                    }
            }

            template <typename site_like, typename mutation_record_like>
            void
            record_site_during_rebuild(const site_like& s, mutation_record_like& mr)
            {
                if (sites.empty() || sites.back().position != s.position)
                    {
                        sites.push_back(s);
                    }
                mr.site = sites.size() - 1;
            }

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
            indexed_edge_container input_left;
            /// The output edge vector. "O" in \cite Kelleher2016-cb, page 13
            indexed_edge_container output_right;
            /// This reflects the length of
            /// tables.edges after last simplification.
            /// It can be used to make sure we only sort
            /// newly-added nodes.
            std::ptrdiff_t edge_offset;
            /// A vector of dead/ancient sample nodes
            std::vector<TS_NODE_INT> preserved_nodes;

            explicit table_collection(const double maxpos)
                : L{maxpos}, nodes{}, edges{}, mutations{}, sites{}, input_left{},
                  output_right{}, edge_offset{0}, preserved_nodes{}
            {
                if (maxpos <= 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument("maxpos must be > 0 and finite");
                    }
            }

            table_collection(const TS_NODE_INT num_initial_nodes,
                             const double initial_time, TS_NODE_INT pop,
                             const double maxpos)
                : L{maxpos}, nodes{}, edges{}, mutations{}, sites{}, input_left{},
                  output_right{}, edge_offset{0}, preserved_nodes{}
            {
                if (maxpos <= 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument("maxpos must be > 0 and finite");
                    }
                for (TS_NODE_INT i = 0; i < num_initial_nodes; ++i)
                    {
                        nodes.push_back(node_t{pop, initial_time});
                    }
            }

            void
            sort_edges() noexcept
            /// Sort the edge table.  On PARENT birth times.
            /// The sorting differs from msprime here. The difference
            /// is that we  assume that birth times are recorded forward in
            /// time rather than backwards.
            ///
            /// In between simplifications, edges are appended to the table
            /// starting at edge_offset.  edges prior to that point are,
            /// by definition, sorted b/c they are the output of simplification.
            /// Thus, we only sort the new edges and then perform a left rotate
            /// w.r.to edge_offset to put the nodes in the correct final order.
            /// The rotation is exactly edges.size() swaps.
            {
                std::sort(edges.begin() + edge_offset, edges.end(),
                          [this](const edge_t& a, const edge_t& b) {
                              auto ga = this->nodes[a.parent].time;
                              auto gb = this->nodes[b.parent].time;
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
                          });
                if (edge_offset > 0)
                    {
                        std::rotate(edges.begin(), edges.begin() + edge_offset,
                                    edges.end());
                    }
#ifndef NDEBUG
                // TODO: allow for exceptions
                // rather than assertions.
                if (edge_offset == 0)
                    {
                        assert(edges_are_sorted());
                    }
                else if (preserved_nodes.empty())
                    {
                        assert(edges_are_sorted());
                    }
                else
                    {
                        assert(edges_are_minimally_sorted());
                    }
#endif
            }

            void
            sort_mutations()
            {
                //mutations are sorted by increasing position
                std::sort(mutations.begin(), mutations.end(),
                          [this](const mutation_t& a, const mutation_t& b) {
                              return sites[a.site].position < sites[b.site].position;
                          });
            }

            void
            sort_tables_for_simplification()
            /// Sorts the tables for simplification, which means only
            /// sorting edge and mutation tables, as the site table
            /// will be rebuilt during simplification.
            {
                sort_edges();
                sort_mutations();
            }

            void
            sort_tables()
            /// Sort all tables.  The site table is rebuilt.
            {
                sort_edges();
                sort_mutations_rebuild_site_table();
            }

            bool
            edges_are_sorted() const noexcept
            /// Test that edges are in proper sort order for simiplification
            /// \version 0.8.0 Change from testing minimal requirement to requirement
            /// for simplification.
            {
                return std::is_sorted(edges.begin(), edges.end(),
                                      [this](const edge_t& a, const edge_t& b) {
                                          auto ga = this->nodes[a.parent].time;
                                          auto gb = this->nodes[b.parent].time;
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
                                      });
            }

            bool
            edges_are_minimally_sorted() const noexcept
            /// Test the MINIMAL sorting requirement.
            /// This minimal condition is important,
            /// as ancient sample tracking will not
            /// guarantee a sort order on the parental
            /// node IDs.
            /// \version 0.8.0 Added to library
            {
                return std::is_sorted(edges.begin(), edges.end(),
                                      [this](const edge_t& a, const edge_t& b) {
                                          auto ga = this->nodes[a.parent].time;
                                          auto gb = this->nodes[b.parent].time;
                                          return ga > gb
                                                 && (std::tie(a.child, a.left)
                                                     < std::tie(b.child, b.left));
                                      });
            }

            void
            clear() noexcept
            /// Clears internal vectors.
            {
                nodes.clear();
                edges.clear();
                mutations.clear();
                preserved_nodes.clear();
                sites.clear();
            }

            void
            record_preserved_nodes(const std::vector<TS_NODE_INT>& node_ids)
            /// Take a list of nodes to record as "ancient samples".
            /// Throws an exception if the nodes are already recorded as such.
            {
                for (auto i : node_ids)
                    {
                        if (static_cast<std::size_t>(i) >= nodes.size())
                            {
                                throw std::invalid_argument(
                                    "node id larger than node table size");
                            }
                        if (std::find(preserved_nodes.begin(), preserved_nodes.end(), i)
                            != preserved_nodes.end())
                            {
                                throw std::invalid_argument("node already recorded as a "
                                                            "preserved_node");
                            }
                        preserved_nodes.push_back(i);
                    }
            }

            std::size_t
            push_back_node(double time, std::int32_t pop)
            {
                nodes.push_back(node_t{pop, time});
                return nodes.size() - 1;
            }

            template <typename... args>
            std::size_t
            emplace_back_node(args&&... Args)
            {
                nodes.emplace_back(node_t{std::forward<args>(Args)...});
                return nodes.size() - 1;
            }

            std::size_t
            push_back_edge(double l, double r, TS_NODE_INT parent, TS_NODE_INT child)
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
                input_left.reserve(edges.size());
                output_right.reserve(edges.size());
                input_left.clear();
                output_right.clear();
                for (auto& e : edges)
                    {
                        assert(e.left < e.right);
                        input_left.emplace_back(e.left, -nodes[e.parent].time, e.parent,
                                                e.child);
                        output_right.emplace_back(e.right, nodes[e.parent].time,
                                                  e.parent, e.child);
                    }
                std::sort(input_left.begin(), input_left.end());
                std::sort(output_right.begin(), output_right.end());
            }

            std::size_t
            register_diploid_offspring(
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const std::int32_t population, const double time)
            /// FIXME: having this here violates ORP
            {
                // TODO: carefully document how to index node times.
                auto next_index = emplace_back_node(population, time);
                if (next_index >= std::numeric_limits<TS_NODE_INT>::max())
                    {
                        throw std::invalid_argument("node index too large");
                    }
                split_breakpoints(breakpoints, parents, next_index);
                return next_index;
            }

            void
            rebuild_site_table()
            /// Complete rebuild of the site table.
            {
                auto site_table_copy(sites);
                sites.clear();
                for (auto& mr : mutations)
                    {
                        auto os = mr.site;
                        record_site_during_rebuild(site_table_copy[mr.site], mr);
                        if (site_table_copy[os].position != sites[mr.site].position)
                            {
                                throw tables_error("error rebuilding site table");
                            }
                    }
            }

            void
            sort_mutations_rebuild_site_table()
            /// O(n*log(n)) time plus O(n) additional memory.
            /// This is called by sort_tables, but also should
            /// be called after adding mutations by some manual
            /// means to a table collection.
            {
                sort_mutations();
                rebuild_site_table();
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

            void
            update_offset()
            {
                edge_offset = edges.size();
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
                       && mutations == b.mutations
                       && preserved_nodes == b.preserved_nodes;
            }

            inline bool
            operator!=(const table_collection& b) const
            {
                return !(*this == b);
            }
        };

        template <typename TableCollectionType, typename mcont_t>
        void
        record_mutations_infinite_sites(
            const TS_NODE_INT u, const mcont_t& mutations,
            const std::vector<std::uint32_t>& new_mutation_keys,
            TableCollectionType& tables)
        /// \version Added in 0.8.0
        /// FIXME: move this to a tree seq recording header
        {
            std::int8_t ancestral_state = 0, derived_state = 1;
            for (auto& k : new_mutation_keys)
                {
                    auto site
                        = tables.emplace_back_site(mutations[k].pos, ancestral_state);
                    if (site >= std::numeric_limits<TS_NODE_INT>::max())
                        {
                            throw std::invalid_argument("site index out of range");
                        }
                    tables.push_back_mutation(u, k, site, derived_state,
                                              mutations[k].neutral);
                }
        }
    } // namespace ts
} // namespace fwdpp

#endif
